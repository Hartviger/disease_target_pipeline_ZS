"""
One-file ChEMBL disease-target script.

It does this:

1. Read disease scope CSV files from 00_summon_the_mondos.py.
2. Collect the MONDO and EFO disease IDs.
3. Write a review text file so you can manually check the IDs.
4. Search ChEMBL for drugs linked to those disease IDs.
5. Keep only human protein targets.
6. Save one ChEMBL target CSV file per disease.

Input:
- disease_scopes/disease_scope_*.csv

Output:
- results/01_chembl_targets/disease_identifier_review.txt
- results/01_chembl_targets/chembl_targets_<disease>.csv

Run:  python3 01_hunt_the_chembl_targets.py
      python3 01_hunt_the_chembl_targets.py focal_epilepsy
"""

import re
import sqlite3
import sys
from pathlib import Path

import pandas as pd


# =============================================================================
# USER INPUT
# =============================================================================

CHEMBL_DB = Path(
    "/Users/jakobhartvig/data/chembl/chembl_37/"
    "chembl_37_sqlite/chembl_37.db"
)

HERE = Path(__file__).resolve().parent
DISEASE_SCOPE_DIR = HERE / "disease_scopes"
RESULTS_DIR = HERE / "results"
CHEMBL_TARGETS_DIR = RESULTS_DIR / "01_chembl_targets"
ID_REVIEW_PATH = CHEMBL_TARGETS_DIR / "disease_identifier_review.txt"

# Leave this as "all" to run every disease_scope_*.csv file in the folder.
# You can also set it to a list, for example: ["focal_epilepsy", "melanoma"].
DISEASES_TO_ANALYZE = "all"

# Creates a simple text file first, so you can manually check the collected IDs.
WRITE_ID_REVIEW = True

# Set this to False if you want to stop after creating the review file.
RUN_CHEMBL_ANALYSIS = True

# 4 = approved drugs only. 1 includes drugs from phase 1 and above.
MINIMUM_PHASE = 4

# 9606 = Homo sapiens / human.
HUMAN_TAX_ID = 9606


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================


def available_diseases():
    """Return the disease names represented by scope files in the folder."""

    prefix = "disease_scope_"
    return sorted(
        path.stem.removeprefix(prefix)
        for path in DISEASE_SCOPE_DIR.glob(f"{prefix}*.csv")
    )


def split_identifiers(value):
    """Split a semicolon-separated CSV cell into individual identifiers."""

    if not value:
        return []

    return [item for item in re.split(r"\s*;\s*", value.strip()) if item]


def normalize_ontology_id(identifier):
    """Normalize MONDO/EFO capitalization for matching against ChEMBL."""

    prefix, separator, number = identifier.partition(":")

    if not separator or prefix.upper() not in {"MONDO", "EFO"}:
        return None

    return f"{prefix.upper()}:{number}"


def selected_diseases():
    """Return the diseases requested by the user settings or command line."""

    if len(sys.argv) > 1:
        return sys.argv[1:]

    if DISEASES_TO_ANALYZE == "all":
        diseases = available_diseases()

        if not diseases:
            raise SystemExit(
                f"No disease scope files found in {DISEASE_SCOPE_DIR}\n"
                "Run 00_summon_the_mondos.py first to create at least one scope file."
            )

        return diseases

    return list(DISEASES_TO_ANALYZE)


def sql_placeholders(values):
    """Return one SQLite placeholder for every supplied value."""

    return ", ".join("?" for _ in values)


def join_unique(values):
    """Join unique, non-empty database values for the final table."""

    return "; ".join(
        sorted({str(value) for value in values if pd.notna(value) and value})
    )


# =============================================================================
# LOAD DISEASE IDENTIFIERS
# =============================================================================


def load_disease_scope(disease):
    """Load all MONDO and EFO identifiers for one disease scope."""

    scope_path = DISEASE_SCOPE_DIR / f"disease_scope_{disease}.csv"

    if not scope_path.exists():
        choices = available_diseases()
        choice_text = ", ".join(choices) if choices else "none found"
        raise SystemExit(
            f"No scope file found for {disease!r}: {scope_path}\n"
            f"Available diseases: {choice_text}\n"
            "Run 00_summon_the_mondos.py first to create the missing scope file."
        )

    scope = pd.read_csv(scope_path, dtype=str, keep_default_na=False)
    required_columns = {"disease_id", "mondo_ids", "efo_ids"}
    missing_columns = sorted(required_columns - set(scope.columns))

    if missing_columns:
        raise SystemExit(
            f"{scope_path.name} is missing columns: {', '.join(missing_columns)}"
        )

    ontology_ids = set()

    for column in ("disease_id", "mondo_ids", "efo_ids"):
        for cell in scope[column]:
            for identifier in split_identifiers(cell):
                normalized = normalize_ontology_id(identifier)

                if normalized:
                    ontology_ids.add(normalized)

    if not ontology_ids:
        raise SystemExit(f"No usable disease identifiers found in {scope_path.name}")

    mondo_ids = sorted(
        identifier for identifier in ontology_ids if identifier.startswith("MONDO:")
    )
    efo_ids = sorted(
        identifier for identifier in ontology_ids if identifier.startswith("EFO:")
    )

    return {
        "disease": disease,
        "scope_path": scope_path,
        "scope": scope,
        "ontology_ids": sorted(ontology_ids),
        "mondo_ids": mondo_ids,
        "efo_ids": efo_ids,
    }


# =============================================================================
# WRITE MANUAL REVIEW FILE # FOR TESTING 
# =============================================================================


def write_identifier_review(loaded_scopes):
    """Write a plain-text checkpoint of the IDs collected from every scope."""

    lines = [
        "Disease identifier review",
        "",
        f"Scope folder: {DISEASE_SCOPE_DIR}",
        "",
        (
            "This file is only for manual checking. ChEMBL will use the "
            "MONDO/EFO IDs as disease ontology IDs."
        ),
        "",
    ]

    for loaded_scope in loaded_scopes:
        scope = loaded_scope["scope"]
        root_rows = scope[scope["relationship"] == "root"]
        root_label = (
            root_rows.iloc[0]["label"]
            if not root_rows.empty
            else "root row not found"
        )

        lines.extend(
            [
                "=" * 78,
                f"Disease name: {loaded_scope['disease']}",
                f"Scope file:   {loaded_scope['scope_path'].name}",
                f"Root label:   {root_label}",
                f"Rows in file: {len(scope)}",
                "",
                f"MONDO IDs ({len(loaded_scope['mondo_ids'])}):",
                ", ".join(loaded_scope["mondo_ids"]) or "none",
                "",
                f"EFO IDs ({len(loaded_scope['efo_ids'])}):",
                ", ".join(loaded_scope["efo_ids"]) or "none",
                "",
                "Terms included:",
            ]
        )

        for row in scope.itertuples(index=False):
            lines.append(
                f"- {row.relationship}: {row.disease_id} | {row.label}"
            )

        lines.append("")

    CHEMBL_TARGETS_DIR.mkdir(parents=True, exist_ok=True)
    ID_REVIEW_PATH.write_text("\n".join(lines), encoding="utf-8")
    return ID_REVIEW_PATH


# =============================================================================
# SEARCH CHEMBL
# =============================================================================


def analyze_disease(loaded_scope):
    """Run the ChEMBL query for one loaded disease scope."""

    disease = loaded_scope["disease"]
    scope_path = loaded_scope["scope_path"]
    ontology_ids = loaded_scope["ontology_ids"]
    query_parameters = ontology_ids + [MINIMUM_PHASE, HUMAN_TAX_ID]

    # One row per drug and target protein, with the phase reached for a matching
    # MONDO or EFO disease identifier.
    query = f"""
    SELECT DISTINCT molecule_dictionary.pref_name       AS drug,
                    drug_indication.max_phase_for_ind   AS phase,
                    drug_indication.efo_id              AS ontology_id,
                    target_dictionary.chembl_id         AS chembl_target_id,
                    target_dictionary.pref_name         AS chembl_target,
                    target_dictionary.target_type       AS chembl_target_type,
                    drug_mechanism.action_type          AS action,
                    component_sequences.tax_id          AS target_tax_id,
                    component_sequences.accession       AS uniprot
    FROM drug_indication
    JOIN molecule_dictionary
        ON drug_indication.molregno = molecule_dictionary.molregno
    JOIN drug_mechanism
        ON drug_mechanism.molregno = drug_indication.molregno
    JOIN target_dictionary
        ON drug_mechanism.tid = target_dictionary.tid
    JOIN target_components
        ON target_dictionary.tid = target_components.tid
    JOIN component_sequences
        ON target_components.component_id = component_sequences.component_id
    WHERE drug_indication.efo_id IN ({sql_placeholders(ontology_ids)})
      AND drug_indication.max_phase_for_ind >= ?
      AND component_sequences.accession IS NOT NULL
      AND component_sequences.tax_id = ?
    ORDER BY uniprot, drug;
    """

    connection = sqlite3.connect(CHEMBL_DB.as_uri() + "?mode=ro", uri=True)

    try:
        drug_target_pairs = pd.read_sql_query(
            query,
            connection,
            params=query_parameters,
        )
    finally:
        connection.close()

    component_counts = (
        drug_target_pairs
        .groupby("chembl_target_id")["uniprot"]
        .nunique()
    )
    drug_target_pairs["target_component_count"] = (
        drug_target_pairs["chembl_target_id"].map(component_counts)
    )

    # Collapse to one row per protein component within one ChEMBL target.
    # This keeps family/complex expansions visible instead of pretending every
    # protein row is an independent pharmacological target.
    targets = (
        drug_target_pairs
        .groupby(["uniprot", "chembl_target_id"])
        .agg(
            ontology_ids=("ontology_id", join_unique),
            chembl_target=("chembl_target", join_unique),
            chembl_target_type=("chembl_target_type", join_unique),
            target_component_count=("target_component_count", "max"),
            action=("action", join_unique),
            target_tax_id=("target_tax_id", join_unique),
            max_phase=("phase", "max"),
            approved_drugs=("drug", join_unique),
            n_drugs=("drug", "nunique"),
        )
        .reset_index()
        .sort_values(["uniprot", "chembl_target_id"])
    )

    print(f"Disease: {disease}")
    print(f"Scope:   {scope_path.name}")
    print(f"Terms:   {len(ontology_ids)} MONDO/EFO ontology IDs")
    print(f"Minimum phase: {MINIMUM_PHASE}")
    print(f"Target organism tax ID: {HUMAN_TAX_ID}")
    print(f"Drug-target rows found: {len(drug_target_pairs)}")
    print(f"Unique drugs found:     {drug_target_pairs['drug'].nunique()}")
    print(f"ChEMBL targets found:   {targets['chembl_target_id'].nunique()}")
    print(f"Protein rows saved:     {len(targets)}")

    CHEMBL_TARGETS_DIR.mkdir(parents=True, exist_ok=True)
    output_path = CHEMBL_TARGETS_DIR / f"chembl_targets_{disease}.csv"
    targets.to_csv(output_path, index=False)
    print(f"Saved CSV: {output_path.name}")

    return output_path


# =============================================================================
# RUN EVERYTHING
# =============================================================================


def main():
    diseases = selected_diseases()
    loaded_scopes = [load_disease_scope(disease) for disease in diseases]

    if WRITE_ID_REVIEW:
        review_path = write_identifier_review(loaded_scopes)
        print(f"Wrote identifier review: {review_path}")

    if not RUN_CHEMBL_ANALYSIS:
        print("RUN_CHEMBL_ANALYSIS is False, so ChEMBL was not searched.")
        return

    for loaded_scope in loaded_scopes:
        print("\n" + "=" * 78)
        analyze_disease(loaded_scope)


if __name__ == "__main__":
    main()
