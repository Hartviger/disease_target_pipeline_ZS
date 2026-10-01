# -----------------------------------------------------------------------------
# Compare ChEMBL approved-drug targets with Open Targets GWAS and RNA scores
# (one row per ChEMBL target). IntOGen scores are saved for stage 04.
# -----------------------------------------------------------------------------

from pathlib import Path
import duckdb


results_folder = Path(__file__).resolve().parent / "results"
open_targets_file = results_folder / "02_open_targets_scores.csv"
chembl_file = Path(
    "/Users/jakobhartvig/ChEMBL_pipeline/results/"
    "03_chembl_target_summary/Unique ChEMBL targets.csv"
)

con = duckdb.connect()

# Split each ChEMBL target into its UniProt components, match them to
# Open Targets, and keep the best score per source for each target.
con.execute(
    f"""
    CREATE TABLE comparison AS
    WITH open_targets AS (
        SELECT
            disease_key,
            symbol,
            UNNEST(STRING_SPLIT(uniprot_ids, ';')) AS uniprot_id,
            gwas_score,
            rna_score,
            intogen_score
        FROM read_csv('{open_targets_file}')
    ),
    chembl AS (
        SELECT
            Disease AS disease,
            "ChEMBL target ID" AS chembl_target_id,
            "ChEMBL target" AS target_name,
            "Target type" AS target_type,
            Gene AS genes,
            UNNEST(STRING_SPLIT(REPLACE("UniProt components", ' ', ''), ';'))
                AS uniprot_id
        FROM read_csv('{chembl_file}')
    )
    SELECT
        c.disease,
        c.chembl_target_id,
        c.target_name,
        c.target_type,
        c.genes,
        COUNT(DISTINCT c.uniprot_id) AS components,
        MAX(o.gwas_score) AS gwas_score,
        ARG_MAX(o.symbol, o.gwas_score) AS gwas_gene,
        MAX(o.rna_score) AS rna_score,
        ARG_MAX(o.symbol, o.rna_score) AS rna_gene,
        MAX(o.intogen_score) AS intogen_score,
        ARG_MAX(o.symbol, o.intogen_score) AS intogen_gene
    FROM chembl c
    LEFT JOIN open_targets o
        ON o.disease_key = REPLACE(LOWER(c.disease), ' ', '_')
       AND o.uniprot_id = c.uniprot_id
    GROUP BY ALL
    """
)

print("ChEMBL targets:", con.sql("SELECT COUNT(*) FROM comparison").fetchone()[0])

# One table per disease (GWAS and RNA only; IntOGen is shown in stage 04).
diseases = con.sql("SELECT DISTINCT disease FROM comparison ORDER BY disease").fetchall()

for (disease,) in diseases:
    print(f"\n{disease}: ChEMBL targets compared with Open Targets")
    con.sql(
        f"""
        SELECT
            chembl_target_id,
            target_name,
            genes,
            ROUND(gwas_score, 4) AS gwas_score,
            gwas_gene,
            ROUND(rna_score, 4) AS rna_score,
            rna_gene
        FROM comparison
        WHERE disease = '{disease}'
        ORDER BY chembl_target_id
        """
    ).show(max_rows=100, max_width=250, max_col_width=40)

# Save everything in one file (scores unrounded, full gene lists).
output_file = results_folder / "03_chembl_target_comparison.csv"
con.execute(f"COPY comparison TO '{output_file}' (HEADER)")
print("\nSaved:", output_file)
