"""
One-file ChEMBL-to-network comparison script.

It does this:

1. Read ChEMBL target CSV files from 01_hunt_the_chembl_targets.py.
2. Convert UniProt IDs into STRING protein IDs.
3. Read the selected disease network from the Cytoscape file.
4. Check if each ChEMBL target maps to STRING.
5. Check if each mapped target is inside the disease network.
6. Add the network level.
7. Save one network comparison CSV file per disease.

Input:
- results/01_chembl_targets/chembl_targets_<disease>.csv
- data/enrichment_and_clustering.cys
- data/9606.protein.aliases.v12.5.txt.gz

Output:
- results/02_network_comparisons/network_comparison_<disease>.csv

Run:  python3 02_map_targets_to_the_network.py
      python3 02_map_targets_to_the_network.py focal_epilepsy
"""

import gzip
import sys
import zipfile
from pathlib import Path

import pandas as pd


# =============================================================================
# USER INPUT
# =============================================================================

HERE = Path(__file__).resolve().parent
RESULTS_DIR = HERE / "results"
CHEMBL_TARGETS_DIR = RESULTS_DIR / "01_chembl_targets"
NETWORK_COMPARISON_DIR = RESULTS_DIR / "02_network_comparisons"
CYS_FILE = HERE / "data" / "enrichment_and_clustering.cys"
ALIASES_FILE = HERE / "data" / "9606.protein.aliases.v12.5.txt.gz"

# Leave this as "all" to compare every chembl_targets_*.csv file from step 01.
# You can also set it to a list, for example: ["focal_epilepsy", "melanoma"].
DISEASES_TO_COMPARE = "all"

# The name each disease network has inside the Cytoscape session.
NETWORK_NAMES = {
    "atopic_dermatitis":  "atopic+dermatitis",
    "ulcerative_colitis": "ulcerative+colitis",
    "focal_epilepsy":     "focal+epilepsy",
    "colorectal_cancer":  "colorectal+cancer",
    "lymphoma":           "lymphoma",
    "melanoma":           "melanoma",
    "aortic_aneurysm":    "aortic+aneurysm",
}


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================


def available_chembl_target_diseases():
    """Return disease names represented by outputs from 01_hunt_the_chembl_targets.py."""

    prefix = "chembl_targets_"
    return sorted(
        path.stem.removeprefix(prefix)
        for path in CHEMBL_TARGETS_DIR.glob(f"{prefix}*.csv")
    )


def selected_diseases():
    """Return the diseases requested by the user settings or command line."""

    if len(sys.argv) > 1:
        diseases = sys.argv[1:]
    elif DISEASES_TO_COMPARE == "all":
        diseases = [
            disease
            for disease in available_chembl_target_diseases()
            if disease in NETWORK_NAMES
        ]
    else:
        diseases = list(DISEASES_TO_COMPARE)

    unknown = sorted(set(diseases) - set(NETWORK_NAMES))

    if unknown:
        sys.exit(
            f"Unknown disease: {', '.join(unknown)}\n"
            f"Choose one of: {', '.join(NETWORK_NAMES)}"
        )

    if not diseases:
        sys.exit(
            f"No chembl_targets_*.csv files found in {CHEMBL_TARGETS_DIR}. "
            "Run 01_hunt_the_chembl_targets.py first."
        )

    return diseases


# =============================================================================
# LOAD DISEASE NETWORK
# =============================================================================


def load_network(disease_name):
    """Read one disease network's node table out of the Cytoscape session."""

    network_name = NETWORK_NAMES[disease_name]

    with zipfile.ZipFile(CYS_FILE) as session:
        # Main networks are stored as "<id>-<name>/...CyNode-<id>+shared+node.cytable".
        # The pathway sub-networks have the name (not the id) in the filename, so
        # requiring a numeric id here keeps us to the full disease networks.
        matches = [
            name for name in session.namelist()
            if f"-{network_name}/SHARED_ATTRS" in name
            and "CyNode" in name
            and name.rsplit("CyNode-", 1)[1].split("+")[0].isdigit()
        ]
        if not matches:
            sys.exit(f"No network table found for {disease_name}")

        with session.open(matches[0]) as table_file:
            # The first rows are Cytoscape metadata, not proteins.
            nodes = pd.read_csv(table_file, skiprows=1).iloc[3:]

    # Remove Cytoscape bookkeeping columns and keep the network annotation we need.
    nodes = nodes.rename(columns={
        "target::development level": "network_level",
    })

    # Prepare the STRING ID for joining with the ChEMBL targets.
    nodes["string_id"] = nodes["@id"].str.replace("stringdb:", "", regex=False)

    return nodes[["string_id", "network_level"]]


# =============================================================================
# LOAD STRING MAPPINGS
# =============================================================================


def load_uniprot_to_string():
    """UniProt accession -> STRING ID, straight from STRING's own alias file."""

    mapping = {}

    with gzip.open(ALIASES_FILE, "rt") as aliases:
        next(aliases)  # header

        for line in aliases:
            string_id, alias, source = line.rstrip("\n").split("\t")

            if source == "UniProt_AC":
                mapping[alias] = string_id

    return mapping


# =============================================================================
# COMPARE CHEMBL TARGETS WITH NETWORK
# =============================================================================


def compare_disease(disease, uniprot_to_string):
    """Build and save one disease-vs-network comparison table."""

    targets = pd.read_csv(CHEMBL_TARGETS_DIR / f"chembl_targets_{disease}.csv")
    required_columns = {
        "uniprot",
        "ontology_ids",
        "chembl_target_id",
        "chembl_target",
        "chembl_target_type",
        "target_component_count",
        "action",
        "max_phase",
        "n_drugs",
        "approved_drugs",
    }
    missing_columns = sorted(required_columns - set(targets.columns))

    if missing_columns:
        sys.exit(
            f"chembl_targets_{disease}.csv is missing columns: "
            f"{', '.join(missing_columns)}\n"
            "Rerun 01_hunt_the_chembl_targets.py before running step 02."
        )

    network = load_network(disease).set_index("string_id")

    targets["string_id"] = targets["uniprot"].map(uniprot_to_string)
    mapped_to_string = targets["string_id"].notna()
    in_network = mapped_to_string & targets["string_id"].isin(network.index)

    network_status = pd.Series("not mapped to STRING", index=targets.index)
    network_status.loc[mapped_to_string] = "not in network"
    network_status.loc[in_network] = "in network"

    comparison = pd.DataFrame({
        "uniprot": targets["uniprot"],
        "ontology_ids": targets["ontology_ids"],
        "chembl_target_id": targets["chembl_target_id"],
        "chembl_target": targets["chembl_target"],
        "chembl_target_type": targets["chembl_target_type"],
        "target_component_count": targets["target_component_count"],
        "action": targets["action"],
        "max_phase": targets["max_phase"],
        "n_drugs": targets["n_drugs"],
        "network_status": network_status,
        "Network level": targets["string_id"].map(network["network_level"]),
        "Approved drugs": targets["approved_drugs"],
    })

    # Keep these as separate states:
    # - not mapped to STRING: we cannot compare it with the network
    # - not in network: it mapped to STRING, but is absent from this disease network
    comparison.loc[~mapped_to_string.values, "Network level"] = "not mapped to STRING"
    comparison.loc[mapped_to_string.values & ~in_network.values, "Network level"] = "not in network"

    status_order = {
        "in network": 0,
        "not in network": 1,
        "not mapped to STRING": 2,
    }
    comparison["_status_order"] = comparison["network_status"].map(status_order)
    comparison = comparison.sort_values(
        ["_status_order", "Network level", "uniprot"],
        ascending=[True, True, True],
    ).drop(columns=["network_status", "_status_order"])

    print(
        f"{disease}: {len(comparison)} targets, "
        f"{targets['chembl_target_id'].nunique()} ChEMBL target IDs, "
        f"{in_network.sum()} in the network, "
        f"{(mapped_to_string & ~in_network).sum()} not in the network, "
        f"{(~mapped_to_string).sum()} not mapped to STRING"
    )

    NETWORK_COMPARISON_DIR.mkdir(parents=True, exist_ok=True)
    output_path = NETWORK_COMPARISON_DIR / f"network_comparison_{disease}.csv"
    comparison.to_csv(output_path, index=False, sep=";", decimal=",")
    print(f"Saved {output_path.name}")

    return output_path


# =============================================================================
# RUN EVERYTHING
# =============================================================================


def main():
    uniprot_to_string = load_uniprot_to_string()

    for disease in selected_diseases():
        compare_disease(disease, uniprot_to_string)


if __name__ == "__main__":
    main()
