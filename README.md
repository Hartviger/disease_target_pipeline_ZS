# disease_target_pipeline_ZS

A co-lab between me, Jakob, student at SDU, and ZS Associates.

## What it does

Finds the approved-drug targets for a disease and checks which of them sit inside that
disease's STRING protein network.

## Pipeline

Run in order. Each stage reads the previous stage's output.

| `00_summon_the_mondos.py` | Resolves a disease name to a MONDO/EFO term via OLS and collects its descendant terms | `disease_scopes/` |

| `01_hunt_the_chembl_targets.py` | Finds approved drugs indicated for those terms, and their human protein targets | `results/01_chembl_targets/` |

| `02_map_targets_to_the_network.py` | Maps UniProt → STRING and checks membership in the disease network | `results/02_network_comparisons/` |

## Setup

```bash
pip install pandas
```


## !!!!
Also needed: **ChEMBL 37 as SQLite**. Stage 01 has the path hardcoded near the top —
edit it to match your machine:
## !!!!

```python
CHEMBL_DB = Path(".../chembl_37_sqlite/chembl_37.db")
```

## Data files

`data/` is excluded from the repo (~58 MB). Recreate it before running stage 02:

| File | Where from |
|---|---|
| `9606.protein.aliases.v12.5.txt.gz` | [STRING downloads](https://string-db.org/cgi/download), taxon 9606 |
| `9606.protein.info.v12.5.txt.gz` | same |
| `enrichment_and_clustering.cys` | Local Cytoscape session, https://zenodo.org/records/19134020 |

Network names inside the `.cys` must match `NETWORK_NAMES` in stage 02.

## Running

```bash
python3 00_summon_the_mondos.py                          # diseases listed inside the script
python3 01_hunt_the_chembl_targets.py                    # every scope file
python3 02_map_targets_to_the_network.py focal_epilepsy   # or just one disease
```
