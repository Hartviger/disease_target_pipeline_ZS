from pathlib import Path
import csv
import requests

import json
from datetime import datetime, timezone


# The folder containing this Python script.
project_folder = Path(__file__).resolve().parent

# The disease scopes were created by the ChEMBL pipeline.
disease_scope_folder = Path(
    "/Users/jakobhartvig/ChEMBL_pipeline/disease_scopes"
)

# Find every disease-scope CSV in the folder.
disease_scope_files = sorted(
    disease_scope_folder.glob("disease_scope_*.csv")
)

# Stop if the folder does not contain any matching files.
if not disease_scope_files:
    raise FileNotFoundError(
        f"No disease-scope files found in {disease_scope_folder}"
    )


print("Number of disease-scope files:", len(disease_scope_files))
print("\nDisease-scope files:")

for file_path in disease_scope_files:
    print("-", file_path.name)



    # Store the root disease from each scope file.
root_diseases = []

for file_path in disease_scope_files:

    # Read the CSV as rows where column names act like dictionary keys.
    with file_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    # Find rows marked as the root disease.
    root_rows = [
        row
        for row in rows
        if row["relationship"].strip().lower() == "root"
    ]

    # Every disease-scope file should contain exactly one root.
    if len(root_rows) != 1:
        raise ValueError(
            f"{file_path.name} contains {len(root_rows)} root rows. "
            "Expected exactly one."
        )

    root_row = root_rows[0]

    # Create a short internal name from the filename.
    disease_key = file_path.stem.removeprefix(
        "disease_scope_"
    )

    # Open Targets uses MONDO_0004980 rather than MONDO:0004980.
    open_targets_id = root_row["disease_id"].strip().replace(
        ":",
        "_",
    )

    root_diseases.append(
        {
            "disease_key": disease_key,
            "disease_name": root_row["label"].strip(),
            "open_targets_id": open_targets_id,
            "scope_file": file_path.name,
        }
    )


print("\nRoot diseases:")

for disease in root_diseases:
    print(
        disease["disease_key"],
        "->",
        disease["open_targets_id"],
        "->",
        disease["disease_name"],
    )



# The public Open Targets GraphQL API.
# The important new part is: associatedTargets(enableIndirect: true)
# This tells Open Targets to count associations for the root disease while applying its rules for appropriate descendant evidence.
api_url = "https://api.platform.opentargets.org/api/v4/graphql"


# Ask Open Targets for the ID and name of one disease.
disease_query = """
query CheckDisease($diseaseId: String!) {
    disease(efoId: $diseaseId) {
        id
        name
        therapeuticAreas {
            id
            name
        }
        associatedTargets(
            enableIndirect: true
            page: {
                index: 0
                size: 1
            }
        ) {
            count
        }
    }
}
"""

# Open Targets therapeutic area "cancer or benign tumor".
cancer_area_id = "MONDO_0045024"



metadata_query = """
query PlatformMetadata {
    meta {
        apiVersion {
            x
            y
            z
        }
        dataVersion {
            year
            month
            iteration
        }
    }
}
"""

def query_open_targets(query, variables):
    """Send one query to the Open Targets API."""

    response = requests.post(
        api_url,
        json={
            "query": query,
            "variables": variables,
        },
        timeout=60,
    )

    # Stop if the web request was unsuccessful.
    response.raise_for_status()

    result = response.json()

    # A GraphQL request can succeed technically but contain API errors.
    if result.get("errors"):
        raise RuntimeError(result["errors"])

    return result["data"]


# Record which Open Targets release is being checked.
metadata_result = query_open_targets(
    metadata_query,
    {},
)

metadata = metadata_result["meta"]

data_version = (
    f"{metadata['dataVersion']['year']}."
    f"{metadata['dataVersion']['month']}"
)

# Some releases have an additional iteration number.
if metadata["dataVersion"]["iteration"] is not None:
    data_version += (
        f".{metadata['dataVersion']['iteration']}"
    )

api_version = (
    f"{metadata['apiVersion']['x']}."
    f"{metadata['apiVersion']['y']}."
    f"{metadata['apiVersion']['z']}"
)

# Check every root disease against Open Targets.
checked_diseases = []

print("\nChecking diseases in Open Targets:")

for disease in root_diseases:

    result = query_open_targets(
        disease_query,
        {
            "diseaseId": disease["open_targets_id"],
        },
    )

    api_disease = result["disease"]

    # A missing result means Open Targets did not recognise the ID.
    if api_disease is None:
        raise ValueError(
            "Open Targets did not recognise "
            f"{disease['open_targets_id']}"
        )

    # Confirm that Open Targets returned the requested identifier.
    if api_disease["id"] != disease["open_targets_id"]:
        raise ValueError(
            f"Requested {disease['open_targets_id']}, "
            f"but Open Targets returned {api_disease['id']}"
        )

    checked_diseases.append(
        {
            **disease,
            "open_targets_name": api_disease["name"],
            "associated_target_count": api_disease[
                "associatedTargets"
            ]["count"],
            "therapeutic_areas": "; ".join(
                area["name"] for area in api_disease["therapeuticAreas"]
            ),
            "is_cancer": any(
                area["id"] == cancer_area_id
                for area in api_disease["therapeuticAreas"]
            ),
        }
    )

    print(
        disease["open_targets_id"],
        "->",
        api_disease["name"],
        "->",
        api_disease["associatedTargets"]["count"],
        "associated targets",
        "->",
        api_disease["therapeuticAreas"],
    )


print(
    f"\nVerified {len(checked_diseases)} "
    "diseases in Open Targets."
)

print("\nOpen Targets data version:", data_version)
print("Open Targets API version:", api_version)

# Create a folder for the pipeline results.
results_folder = project_folder / "results"
results_folder.mkdir(parents=True, exist_ok=True)

# Save the verified disease mappings for the download script.
mapping_file = results_folder / "00_disease_mapping.csv"

with mapping_file.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=list(checked_diseases[0].keys()),
    )

    writer.writeheader()
    writer.writerows(checked_diseases)

print("\nDisease mapping saved to:", mapping_file)


# Describe when and how we checked Open Targets.
run_metadata = {
    "checked_at_utc": datetime.now(timezone.utc).isoformat(),
    "api_url": api_url,
    "data_version": data_version,
    "api_version": api_version,
    "number_of_diseases": len(checked_diseases),
    "enable_indirect": True,
    "count_scope": "all evidence sources",
}

# Save this information beside the disease mapping.
metadata_file = results_folder / "00_run_metadata.json"

with metadata_file.open("w", encoding="utf-8") as handle:
    json.dump(run_metadata, handle, indent=4)
    handle.write("\n")

print("Run metadata saved to:", metadata_file)

print("Number of disease-scope files:", len(disease_scope_files))
print("\nDisease-scope files:")

for file_path in disease_scope_files:
    print("-", file_path.name)


    print("\nRoot diseases:")

for disease in root_diseases:
    print(
        disease["disease_key"],
        "->",
        disease["open_targets_id"],
        "->",
        disease["disease_name"],
    )