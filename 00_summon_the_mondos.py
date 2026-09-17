"""
One-file OLS disease scope script.

It does this:
1. Take a disease name or ID.
2. Search OLS.
3. Pick one disease term.
4. Get its descendants.
5. Extract MONDO and EFO IDs.
6. Save a CSV file.
"""

import csv
import json
import re
from pathlib import Path
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen


# =============================================================================
# USER INPUT
# =============================================================================

# Add as many diseases as you want here.
# Each disease gets its own CSV file.
DISEASES = [
    "atopic dermatitis",
    "aortic aneurysm",
    "colorectal cancer",
    "focal epilepsy",
    "lymphoma",
    "melanoma",
    "ulcerative colitis",

    "breast cancer",
]
# example: 
#    "focal epilepsy",
#    "atopic dermatitis",
#    "MONDO:0006526",

# the search will be able to undertand all these search 
# Exact Synonyms: AD / allergic form of dermatitis / ATOD / atopic dermatitis / atopic eczema / dermatitis, atopic / eczema, atopic / eczema / eczematous dermatitis 
# Beware that other diseases have similar synonyms such as "AD"
# MONDO is preferred, because this pipeline starts from disease ontology terms.
# EFO is kept as a fallback if MONDO does not find the disease name.

SEARCH_ONTOLOGIES = ("mondo", "efo")

# Output folder used by the next scripts in the pipeline.
OUTPUT_FOLDER = Path(__file__).resolve().parent / "disease_scopes"


# =============================================================================
# OLS API HELPER
# =============================================================================

OLS_BASE = "https://www.ebi.ac.uk/ols4/api"


def get_json(url, params=None):
    """Ask OLS for data and return it as Python dictionaries/lists."""

    if params:
        url = f"{url}?{urlencode(params)}"

    request = Request(url, headers={"User-Agent": "simple-ols-scope/0.1"})

    with urlopen(request, timeout=60) as response:
        return json.load(response)


def curie_to_iri(curie):
    """Convert MONDO:0005384 or EFO:0004263 into the full OLS IRI."""

    prefix, number = curie.upper().split(":", 1)

    if prefix == "MONDO":
        return "http://purl.obolibrary.org/obo/MONDO_" + number

    if prefix == "EFO":
        return "http://www.ebi.ac.uk/efo/EFO_" + number

    raise ValueError(f"Only MONDO and EFO IDs are supported here: {curie}")


def get_term(ontology, iri):
    """Get one exact OLS term from its ontology and IRI."""

    data = get_json(
        f"{OLS_BASE}/ontologies/{ontology}/terms",
        {"iri": iri},
    )
    terms = data.get("_embedded", {}).get("terms", [])

    if not terms:
        raise LookupError(f"Could not find {iri} in {ontology.upper()}")

    return terms[0]


# =============================================================================
# SEARCH OLS AND PICK DISEASE
# =============================================================================


def normal_text(text):
    """Make text easier to compare."""

    return " ".join(str(text).casefold().split())


def search_ols(query):
    """Search OLS by disease name."""

    data = get_json(
        f"{OLS_BASE}/search",
        {
            "q": query,
            "ontology": ",".join(SEARCH_ONTOLOGIES),
            "type": "class",
            "queryFields": "label,synonym",
            "fieldList": "iri,label,obo_id,ontology_name,synonym",
            "exact": "true",
            "local": "true",
            "rows": 25,
        },
    )

    return data.get("response", {}).get("docs", [])


def pick_disease(query):
    """Pick the best OLS disease result."""

    # If the user typed a real ontology ID, use that directly.
    if re.fullmatch(r"(?i)(MONDO|EFO):\d+", query.strip()):
        ontology = query.split(":", 1)[0].lower()
        term = get_term(ontology, curie_to_iri(query))
        return term, ontology, "ID search"

    results = search_ols(query)

    if not results:
        raise LookupError(f'No exact OLS match found for "{query}".')

    ontology_rank = {name: number for number, name in enumerate(SEARCH_ONTOLOGIES)}
    wanted = normal_text(query)

    def score(result):
        ontology = result.get("ontology_name", "")
        label = result.get("label", "")
        exact_label = normal_text(label) == wanted

        # Lower score is better:
        # 1. Prefer MONDO over EFO.
        # 2. Prefer exact label over synonym match.
        return (
            ontology_rank.get(ontology, 99),
            0 if exact_label else 1,
        )

    selected = sorted(results, key=score)[0]
    ontology = selected["ontology_name"]
    term = get_term(ontology, selected["iri"])
    return term, ontology, "OLS name search"


# =============================================================================
# GET DESCENDANTS
# =============================================================================


def get_descendants(ontology, iri):
    """Get all child disease terms below the selected disease."""

    # OLS needs the IRI to be encoded twice in this endpoint.
    encoded_iri = quote(quote(iri, safe=""), safe="")
    url = f"{OLS_BASE}/ontologies/{ontology}/terms/{encoded_iri}/hierarchicalDescendants"

    descendants = []
    page = 0

    while True:
        data = get_json(url, {"page": page, "size": 500})  # OLS allows max 500 results per page.
        descendants.extend(data.get("_embedded", {}).get("terms", []))

        total_pages = data.get("page", {}).get("totalPages", 1)
        page += 1  # Move to the next page if there are more descendants.

        if page >= total_pages:
            break

    return descendants


# =============================================================================
# EXTRACT MONDO AND EFO IDS
# =============================================================================


def listify(value):
    """OLS sometimes gives one value, sometimes a list. This makes it a list."""

    if value is None:
        return []

    if isinstance(value, list):
        return value

    return [value]


def extract_ids(term, relationship):
    """Turn one OLS term into one CSV row."""

    disease_id = term.get("obo_id", "")

    xrefs = listify(
        term.get("annotation", {}).get("database_cross_reference")
    )

    all_ids = {disease_id, *xrefs}

    mondo_ids = sorted(x for x in all_ids if str(x).startswith("MONDO:"))
    efo_ids = sorted(x for x in all_ids if str(x).startswith("EFO:"))
    return {
        "relationship": relationship,
        "disease_id": disease_id,
        "label": term.get("label", ""),
        "mondo_ids": ";".join(mondo_ids),
        "efo_ids": ";".join(efo_ids),
    }


def build_rows(root_term, descendants):
    """Combine the root disease and all descendants into CSV rows."""

    rows = [extract_ids(root_term, "root")]
    rows.extend(extract_ids(term, "descendant") for term in descendants)

    unique_rows = []
    seen_disease_ids = set()

    for row in rows:
        disease_id = row["disease_id"]

        if disease_id and disease_id not in seen_disease_ids:
            unique_rows.append(row)
            seen_disease_ids.add(disease_id)

    return unique_rows


# =============================================================================
# SAVE CSV
# =============================================================================


def safe_filename(name):
    """Make a safe filename from the disease name."""

    slug = re.sub(r"[^a-z0-9]+", "_", name.casefold()).strip("_")
    return f"disease_scope_{slug or 'disease'}.csv"


def save_csv(rows, disease_name):
    """Save the disease scope table."""

    OUTPUT_FOLDER.mkdir(exist_ok=True)
    output_path = OUTPUT_FOLDER / safe_filename(disease_name)

    columns = [
        "relationship",
        "disease_id",
        "label",
        "mondo_ids",
        "efo_ids",
    ]

    with output_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)

    return output_path


# =============================================================================
# RUN EVERYTHING
# =============================================================================


def run_one_disease(disease_query):
    """Run the full OLS scope process for one disease."""

    print("")
    print(f"Search OLS: {disease_query}")

    selected_term, ontology, method = pick_disease(disease_query)

    print(
        "Selected:",
        selected_term.get("obo_id"),
        "-",
        selected_term.get("label"),
        f"({ontology.upper()}, {method})",
    )

    print("Getting descendants...")
    descendants = get_descendants(ontology, selected_term["iri"])

    print("Extracting MONDO/EFO IDs...")
    rows = build_rows(selected_term, descendants)

    output_path = save_csv(rows, disease_query)

    print(f"Saved {len(rows)} rows to:")
    print(output_path)


def main():
    for disease_query in DISEASES:
        try:
            run_one_disease(disease_query)
        except Exception as error:
            print("")
            print(f"Could not process {disease_query}: {error}")


if __name__ == "__main__":
    main()
