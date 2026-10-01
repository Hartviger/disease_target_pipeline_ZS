# -----------------------------------------------------------------------------
# Sensitivity analysis (cancers only): does somatic evidence from IntOGen
# explain ChEMBL targets that GWAS misses? Reuses the stage 03 output.
# -----------------------------------------------------------------------------

from pathlib import Path
import duckdb


results_folder = Path(__file__).resolve().parent / "results"
comparison_file = results_folder / "03_chembl_target_comparison.csv"
mapping_file = results_folder / "00_disease_mapping.csv"

con = duckdb.connect()

# Keep only the diseases Open Targets classifies as cancer (stage 00).
con.execute(
    f"""
    CREATE TABLE cancer_comparison AS
    SELECT c.*
    FROM read_csv('{comparison_file}') c
    JOIN read_csv('{mapping_file}') m
        ON m.disease_key = REPLACE(LOWER(c.disease), ' ', '_')
    WHERE CAST(m.is_cancer AS BOOLEAN)
    """
)

cancers = con.sql(
    "SELECT DISTINCT disease FROM cancer_comparison ORDER BY disease"
).fetchall()

if not cancers:
    raise ValueError("No cancers found. Check is_cancer in 00_disease_mapping.csv.")

# One table per cancer.
for (disease,) in cancers:
    print(f"\n{disease}: ChEMBL targets, GWAS vs IntOGen")
    con.sql(
        f"""
        SELECT
            chembl_target_id,
            target_name,
            genes,
            ROUND(gwas_score, 4) AS gwas_score,
            gwas_gene,
            ROUND(intogen_score, 4) AS intogen_score,
            intogen_gene
        FROM cancer_comparison
        WHERE disease = '{disease}'
        ORDER BY chembl_target_id
        """
    ).show(max_rows=100, max_width=250, max_col_width=40)

output_file = results_folder / "04_intogen_cancer_sensitivity.csv"
con.execute(f"COPY cancer_comparison TO '{output_file}' (HEADER)")
print("\nSaved:", output_file)
