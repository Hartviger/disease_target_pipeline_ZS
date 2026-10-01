# -----------------------------------------------------------------------------
# Stage 01: Download the Open Targets Parquet files for the release
# recorded in stage 00.
#
# Files are saved in data/<release>/<dataset>/...
# A manifest with file sizes and SHA-256 checksums is saved in
# results/01_download_manifest.json.
# Files that are already complete are not downloaded again.
# -----------------------------------------------------------------------------

from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
import shutil
import requests


# -----------------------------------------------------------------------------
# Setup: paths and release from stage 00
# -----------------------------------------------------------------------------

project_folder = Path(__file__).resolve().parent
results_folder = project_folder / "results"
metadata_file = results_folder / "00_run_metadata.json"
manifest_file = results_folder / "01_download_manifest.json"

if not metadata_file.is_file():
    raise FileNotFoundError(
        f"Required input is missing: {metadata_file}. "
        "Run 00_check_open_targets.py first."
    )

with metadata_file.open(encoding="utf-8") as handle:
    check_metadata = json.load(handle)

data_version = check_metadata["data_version"]

release_url = (
    "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/"
    f"{data_version}/output/"
)

data_folder = project_folder / "data" / data_version

# Datasets needed by the pipeline:
# - association_overall_indirect: overall score per disease-target pair
# - association_by_datasource_indirect: score per datasource
#   (gwas_credible_sets, expression_atlas, ...)
# - target: gene annotation, used to map Ensembl IDs to UniProt
datasets = [
    "association_overall_indirect",
    "association_by_datasource_indirect",
    "target",
]

session = requests.Session()

print("Release:", data_version)
print("Server folder:", release_url)
print("Local folder:", data_folder)


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

def list_links(url):
    """Return the file and folder names shown in a server directory listing."""

    response = session.get(url, timeout=60)
    response.raise_for_status()

    links = re.findall(r'href="([^"?]+)"', response.text)

    return [
        link for link in links
        if not link.startswith(("/", "http", ".."))
    ]


def find_parquet_files(url, depth=0):
    """Find all Parquet files in a folder, including subfolders."""

    parquet_files = []

    for link in list_links(url):
        if link.endswith(".parquet"):
            parquet_files.append(url + link)
        elif link.endswith("/") and depth < 3:
            parquet_files.extend(find_parquet_files(url + link, depth + 1))

    return parquet_files


def remote_size(url):
    """Ask the server for the size of one file (bytes) without downloading."""

    response = session.head(url, allow_redirects=True, timeout=60)
    response.raise_for_status()
    size = response.headers.get("Content-Length")

    if size is None:
        raise RuntimeError(f"The server did not report a size for {url}")

    return int(size)


def sha256_of_file(path):
    """Calculate the SHA-256 checksum of a local file."""

    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def download_file(url, destination, expected_size):
    """Download one file. Write to a .part file first, then rename."""

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".part")
    digest = hashlib.sha256()

    with session.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()

        with temporary.open("wb") as handle:
            for block in response.iter_content(chunk_size=1024 * 1024):
                handle.write(block)
                digest.update(block)

    if temporary.stat().st_size != expected_size:
        raise RuntimeError(
            f"Incomplete download: {url} "
            f"({temporary.stat().st_size} of {expected_size} bytes)"
        )

    temporary.replace(destination)
    return digest.hexdigest()


def readable(size_in_bytes):
    """Show a byte count as MB or GB."""

    if size_in_bytes >= 1024 ** 3:
        return f"{size_in_bytes / 1024 ** 3:.2f} GB"
    return f"{size_in_bytes / 1024 ** 2:.1f} MB"


# -----------------------------------------------------------------------------
# Step 1: list the files and their sizes
# -----------------------------------------------------------------------------

output_entries = list_links(release_url)
files_to_get = []

for dataset in datasets:
    if f"{dataset}/" not in output_entries:
        relevant = [
            entry for entry in output_entries
            if "association" in entry or entry.startswith("target")
        ]
        raise ValueError(
            f"Dataset '{dataset}' was not found in release {data_version}. "
            f"Relevant folders on the server: {relevant}"
        )

    dataset_url = release_url + dataset + "/"
    parquet_urls = find_parquet_files(dataset_url)

    if not parquet_urls:
        raise ValueError(f"No Parquet files found for {dataset}")

    dataset_size = 0

    for url in parquet_urls:
        size = remote_size(url)
        dataset_size += size

        files_to_get.append(
            {
                "dataset": dataset,
                "url": url,
                "relative_path": f"{dataset}/{url.removeprefix(dataset_url)}",
                "size_bytes": size,
            }
        )

    print(f"\n{dataset}: {len(parquet_urls)} files, {readable(dataset_size)}")


total_size = sum(item["size_bytes"] for item in files_to_get)
free_space = shutil.disk_usage(project_folder).free

print("\nTotal size:", readable(total_size))
print("Free disk space:", readable(free_space))

# Keep a safety margin of 5 GB.
if total_size + 5 * 1024 ** 3 > free_space:
    raise RuntimeError("Not enough free disk space for the download.")

answer = input("\nDownload these files? [y/N] ").strip().lower()

if answer != "y":
    raise SystemExit("Download cancelled. Nothing was downloaded.")


# -----------------------------------------------------------------------------
# Step 2: download (skip files that are already complete)
# -----------------------------------------------------------------------------

# Reuse checksums from an earlier run so existing files are not re-read.
previous_checksums = {}

if manifest_file.is_file():
    with manifest_file.open(encoding="utf-8") as handle:
        previous_manifest = json.load(handle)

    if previous_manifest.get("data_version") == data_version:
        previous_checksums = {
            item["relative_path"]: item["sha256"]
            for item in previous_manifest["files"]
        }

for number, item in enumerate(files_to_get, start=1):
    destination = data_folder / item["relative_path"]
    label = f"[{number}/{len(files_to_get)}] {item['relative_path']}"

    already_complete = (
        destination.is_file()
        and destination.stat().st_size == item["size_bytes"]
    )

    if already_complete:
        item["sha256"] = (
            previous_checksums.get(item["relative_path"])
            or sha256_of_file(destination)
        )
        print(label, "already present")
        continue

    item["sha256"] = download_file(
        item["url"], destination, item["size_bytes"]
    )
    print(label, readable(item["size_bytes"]))


# -----------------------------------------------------------------------------
# Step 3: save the manifest
# -----------------------------------------------------------------------------

manifest = {
    "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
    "data_version": data_version,
    "release_url": release_url,
    "local_folder": str(data_folder),
    "datasets": datasets,
    "number_of_files": len(files_to_get),
    "total_size_bytes": total_size,
    "files": files_to_get,
}

with manifest_file.open("w", encoding="utf-8") as handle:
    json.dump(manifest, handle, indent=4)
    handle.write("\n")

print("\nAll files downloaded and checked.")
print("Manifest saved to:", manifest_file)


# -----------------------------------------------------------------------------
# Step 4: quick look with DuckDB (columns and row counts)
# -----------------------------------------------------------------------------

try:
    import duckdb
except ImportError:
    print("\nDuckDB is not installed. Run: pip install duckdb")
    raise SystemExit

for dataset in datasets:
    pattern = str(data_folder / dataset / "**" / "*.parquet")

    print(f"\n{dataset}:")
    print(duckdb.sql(f"DESCRIBE SELECT * FROM read_parquet('{pattern}')"))
    print(duckdb.sql(f"SELECT COUNT(*) AS rows FROM read_parquet('{pattern}')"))
