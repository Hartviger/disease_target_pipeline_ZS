
# -----------------------------------------------------------------------------
# Setup: libraries, file paths, and inputs from stages 00 and 01
# -----------------------------------------------------------------------------

from pathlib import Path
import csv
import json
import duckdb


# The folder containing this Python script.
project_folder = Path(__file__).resolve().parent
results_folder = project_folder / "results"

# Read the data release recorded in stage 00.
with (results_folder / "00_run_metadata.json").open(encoding="utf-8") as handle:
    data_version = json.load(handle)["data_version"]

# The Parquet files downloaded in stage 01.
data_folder = project_folder / "data" / data_version
by_datasource_files = (
    data_folder / "association_by_datasource_indirect" / "*.parquet"
)

# Read the seven diseases verified in stage 00.
with (results_folder / "00_disease_mapping.csv").open(
    newline="", encoding="utf-8"
) as handle:
    diseases = list(csv.DictReader(handle))


print("Data release:", data_version)
print("Parquet folder:", data_folder)
print("Number of diseases:", len(diseases))

print("\nDiseases:")

# -----------------------------------------------------------------------------
# Read GWAS and Expression Atlas scores for the seven diseases
# -----------------------------------------------------------------------------

for disease in diseases:
    print("-", disease["disease_key"], "->", disease["open_targets_id"])



# The Open Targets datasources we want to keep (IntOGen is used in stage 04).
datasources = ["gwas_credible_sets", "expression_atlas", "intogen"]

# Turn the Python lists into SQL lists, e.g. 'MONDO_0005160', 'MONDO_0004980'.
disease_id_list = ", ".join(
    f"'{disease['open_targets_id']}'" for disease in diseases
)
datasource_list = ", ".join(f"'{name}'" for name in datasources)

# An in-memory DuckDB database (nothing is written to disk).
con = duckdb.connect()

# Keep only our seven diseases and the two datasources.
con.execute(
    f"""
    CREATE TABLE long_scores AS
    SELECT
        diseaseId,
        targetId,
        aggregationValue AS datasource,
        associationScore AS score
    FROM read_parquet('{by_datasource_files}')
    WHERE aggregationType = 'datasourceId'
      AND diseaseId IN ({disease_id_list})
      AND aggregationValue IN ({datasource_list})
    """
)

# How many targets does each disease have per datasource?
print("\nTargets per disease and datasource:")
print(
    con.sql(
        """
        SELECT diseaseId, datasource, COUNT(*) AS targets
        FROM long_scores
        GROUP BY diseaseId, datasource
        ORDER BY diseaseId, datasource
        """
    )
)

# -----------------------------------------------------------------------------
# Combine into one row per disease-gene pair
# -----------------------------------------------------------------------------

# Safety check: each disease-target-datasource combination should appear once.
duplicate_rows = con.sql(
    """
    SELECT COUNT(*) - COUNT(DISTINCT (diseaseId, targetId, datasource))
    FROM long_scores
    """
).fetchone()[0]

if duplicate_rows:
    raise RuntimeError(f"{duplicate_rows} duplicate rows in long_scores.")


# One row per disease-target pair, one score column per datasource.
con.execute(
    """
    CREATE TABLE pair_scores AS
    SELECT
        diseaseId,
        targetId,
        MAX(score) FILTER (WHERE datasource = 'gwas_credible_sets') AS gwas_score,
        MAX(score) FILTER (WHERE datasource = 'expression_atlas') AS rna_score,
        MAX(score) FILTER (WHERE datasource = 'intogen') AS intogen_score
    FROM long_scores
    GROUP BY diseaseId, targetId
    """
)

# How many genes have GWAS, RNA or both?
print("\nDisease-gene pairs per disease:")
print(
    con.sql(
        """
        SELECT
            diseaseId,
            COUNT(*) AS genes,
            COUNT(gwas_score) AS with_gwas,
            COUNT(rna_score) AS with_rna,
            COUNT(*) FILTER (
                WHERE gwas_score IS NOT NULL AND rna_score IS NOT NULL
            ) AS with_both
        FROM pair_scores
        GROUP BY diseaseId
        ORDER BY diseaseId
        """
    )
)


# -----------------------------------------------------------------------------
# Build the Ensembl-to-UniProt lookup (Swiss-Prot only)
# -----------------------------------------------------------------------------

# The target dataset downloaded in stage 01.
target_files = data_folder / "target" / "*.parquet"

# Which kinds of protein IDs does Open Targets store?
print("\nProtein ID sources in the target dataset:")
print(
    con.sql(
        f"""
        SELECT protein.source AS source, COUNT(*) AS ids
        FROM (
            SELECT UNNEST(proteinIds) AS protein
            FROM read_parquet('{target_files}')
        )
        GROUP BY source
        ORDER BY ids DESC
        """
    )
)

# One row per gene and reviewed (Swiss-Prot) UniProt accession.
# Genes without a Swiss-Prot accession are kept, with uniprot_id empty.
con.execute(
    f"""
    CREATE TABLE target_uniprot AS
    WITH proteins AS (
        SELECT id, UNNEST(proteinIds) AS protein
        FROM read_parquet('{target_files}')
    ),
    swissprot AS (
        SELECT DISTINCT id, protein.id AS uniprot_id
        FROM proteins
        WHERE protein.source = 'uniprot_swissprot'
    )
    SELECT
        t.id AS targetId,
        t.approvedSymbol AS symbol,
        t.biotype,
        s.uniprot_id
    FROM read_parquet('{target_files}') t
    LEFT JOIN swissprot s ON s.id = t.id
    """
)

# Stop if the Swiss-Prot filter found nothing (e.g. a different source name).
swissprot_rows = con.sql(
    "SELECT COUNT(uniprot_id) FROM target_uniprot"
).fetchone()[0]

if swissprot_rows == 0:
    raise RuntimeError(
        "No 'uniprot_swissprot' IDs found. Check the source names above."
    )

print("Genes in target dataset:",
      con.sql("SELECT COUNT(DISTINCT targetId) FROM target_uniprot").fetchone()[0])
print("Gene-to-Swiss-Prot links:", swissprot_rows)


# -----------------------------------------------------------------------------
# Check the UniProt translation for our disease-gene pairs
# -----------------------------------------------------------------------------

# Per gene in our pairs: is it in the target dataset, and how many
# Swiss-Prot accessions does it have?
print("\nUniProt translation of the Open Targets genes in our pairs:")
print(
    con.sql(
        """
        WITH genes AS (
            SELECT DISTINCT targetId FROM pair_scores
        ),
        per_gene AS (
            SELECT
                g.targetId,
                COUNT(tu.targetId) > 0 AS in_target_dataset,
                COUNT(tu.uniprot_id) AS n_uniprot
            FROM genes g
            LEFT JOIN target_uniprot tu USING (targetId)
            GROUP BY g.targetId
        )
        SELECT
            COUNT(*) AS genes,
            COUNT(*) FILTER (WHERE NOT in_target_dataset) AS not_in_target_dataset,
            COUNT(*) FILTER (WHERE n_uniprot = 0) AS no_swissprot,
            COUNT(*) FILTER (WHERE n_uniprot = 1) AS one_swissprot,
            COUNT(*) FILTER (WHERE n_uniprot > 1) AS several_swissprot
        FROM per_gene
        """
    )
)

# What kind of genes have no Swiss-Prot accession?
print("\nBiotypes of genes without a Swiss-Prot accession:")
print(
    con.sql(
        """
        SELECT tu.biotype, COUNT(DISTINCT tu.targetId) AS genes
        FROM target_uniprot tu
        WHERE tu.targetId IN (SELECT targetId FROM pair_scores)
          AND tu.uniprot_id IS NULL
        GROUP BY tu.biotype
        ORDER BY genes DESC
        """
    )
)

# Accessions shared by more than one gene (in the whole target dataset).
shared_accessions = con.sql(
    """
    SELECT COUNT(*)
    FROM (
        SELECT uniprot_id
        FROM target_uniprot
        WHERE uniprot_id IS NOT NULL
        GROUP BY uniprot_id
        HAVING COUNT(DISTINCT targetId) > 1
    )
    """
).fetchone()[0]

print("\nSwiss-Prot accessions linked to more than one gene:", shared_accessions)


# -----------------------------------------------------------------------------
# Build the final table: one row per disease-gene pair, and save it
# -----------------------------------------------------------------------------

# A small table with the disease key and name for each Open Targets ID.
con.execute(
    """
    CREATE TABLE disease_names (
        diseaseId VARCHAR,
        disease_key VARCHAR,
        disease_name VARCHAR
    )
    """
)

con.executemany(
    "INSERT INTO disease_names VALUES (?, ?, ?)",
    [
        (d["open_targets_id"], d["disease_key"], d["open_targets_name"])
        for d in diseases
    ],
)

# Add disease names, gene symbol, biotype and UniProt IDs to every pair.
# Genes with several accessions get them in one cell, separated by ";".
con.execute(
    """
    CREATE TABLE prepared AS
    SELECT
        d.disease_key,
        d.disease_name,
        p.diseaseId,
        p.targetId,
        ANY_VALUE(tu.symbol) AS symbol,
        ANY_VALUE(tu.biotype) AS biotype,
        STRING_AGG(tu.uniprot_id, ';' ORDER BY tu.uniprot_id) AS uniprot_ids,
        p.gwas_score,
        p.rna_score,
        p.intogen_score
    FROM pair_scores p
    JOIN disease_names d USING (diseaseId)
    LEFT JOIN target_uniprot tu USING (targetId)
    GROUP BY
        d.disease_key, d.disease_name, p.diseaseId, p.targetId,
        p.gwas_score, p.rna_score, p.intogen_score
    ORDER BY d.disease_key, p.targetId
    """
)

# The final table must have exactly one row per pair.
pair_count = con.sql("SELECT COUNT(*) FROM pair_scores").fetchone()[0]
prepared_count = con.sql("SELECT COUNT(*) FROM prepared").fetchone()[0]

if prepared_count != pair_count:
    raise RuntimeError(
        f"Expected {pair_count} rows, but the final table has {prepared_count}."
    )

# Save as CSV.
output_file = results_folder / "02_open_targets_scores.csv"
con.execute(f"COPY prepared TO '{output_file}' (HEADER, DELIMITER ',')")

print("\nRows saved:", prepared_count)
print("Output file:", output_file)
print("\nFirst rows:")
print(con.sql("SELECT * FROM prepared LIMIT 5"))