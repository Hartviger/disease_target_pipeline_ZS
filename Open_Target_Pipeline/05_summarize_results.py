# -----------------------------------------------------------------------------
# Summary: ChEMBL targets with their Open Targets GWAS, RNA and IntOGen scores
# -----------------------------------------------------------------------------

from pathlib import Path
import duckdb


results_folder = Path(__file__).resolve().parent / "results"
comparison_file = results_folder / "03_chembl_target_comparison.csv"
chembl_file = Path(
    "/Users/jakobhartvig/ChEMBL_pipeline/results/"
    "03_chembl_target_summary/Unique ChEMBL targets.csv"
)
output_file = results_folder / "05_summary.csv"

con = duckdb.connect()

# Four columns from ChEMBL, plus the best score per source from stage 03.
con.execute(
    f"""
    COPY (
        SELECT
            t.Disease,
            t."ChEMBL target ID",
            t."Target network status",
            ROUND(c.gwas_score, 4) AS "GWAS score",
            ROUND(c.rna_score, 4) AS "RNA score",
            ROUND(c.intogen_score, 4) AS "IntOGen score",
            T.Gene
        FROM read_csv('{chembl_file}') t
        LEFT JOIN read_csv('{comparison_file}') c
            ON c.disease = t.Disease
           AND c.chembl_target_id = t."ChEMBL target ID"
        ORDER BY t.Disease, t."ChEMBL target ID"
    ) TO '{output_file}' (HEADER)
    """
)

rows = con.sql(f"SELECT COUNT(*) FROM read_csv('{output_file}')").fetchone()[0]
print("Rows saved:", rows)
print("Output file:", output_file)