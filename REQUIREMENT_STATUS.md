# Requirement Status

## Current Policy

Status: updated.

Because the `2 positive + 1 negative` version produced fewer than 10k candidate slots, the active ChEMBL-side dataset now uses:

```text
1 input + exactly 3 property objective groups -> exactly 2 positive answers
```

Negative candidates are omitted from the active ChEMBL-side output.

## Counts

| File | Queries | Positive answer slots | Unique inputs | Primary endpoints | Preserved secondary endpoints |
|---|---:|---:|---:|---:|---:|
| `data/chembl_3prop_2pos/all.csv` | 6,070 | 12,140 | 3,994 | 237 | 261 |
| `data/chembl_3prop_2pos/train.csv` | 5,858 | 11,716 | 3,803 | 235 | 258 |
| `data/chembl_3prop_2pos/val.csv` | 104 | 208 | 99 | 23 | 22 |
| `data/chembl_3prop_2pos/test.csv` | 108 | 216 | 92 | 29 | 28 |

The ChEMBL release layout now matches PubChem: `all.csv`, `train.csv`, `val.csv`, and `test.csv`. The former `strict_all.csv` view was renamed to `all.csv`; the older loose `pretrain.csv` view was folded into `train.csv` after filtering out held-out molecule overlap.

## Other Sources

Processed source-audit files are in:

```text
data/remaining_sources_3prop_2pos/
```

PubChem BioAssay is now also normalized into a separate experimental observation folder:

| File | Content |
|---|---|
| `data/pubchem_bioassay/pubchem_bioassay_observations.csv` | 127,360 pilot normalized PubChem molecule-level assay observations from the initial small shard batch. |
| `data/pubchem_bioassay/pubchem_property_summary.csv` | 27,831 PubChem assay-local endpoint summaries. |
| `data/pubchem_bioassay/pubchem_bioassay_stats.json` | PubChem processing stats and counters. |
| `data/pubchem_bioassay/pubchem_final_supported_observations_merged.csv` | 2,879,222 full-crawl exact p-scale PubChem observations from all 1,933 BioAssay CSV/Data shards. |

PubChem source-specific multi-answer rows are stored separately:

| File | Queries | Positive answers |
|---|---:|---:|
| `data/pubchem_3prop_2pos/all.csv` | 42,754 | 85,508 |
| `data/pubchem_3prop_2pos/train.csv` | 34,204 | 68,408 |
| `data/pubchem_3prop_2pos/val.csv` | 4,275 | 8,550 |
| `data/pubchem_3prop_2pos/test.csv` | 4,275 | 8,550 |

The PubChem split has zero molecule overlap between train/val/test according to `data/pubchem_3prop_2pos/dataset_stats.json`.

Papyrus source-specific multi-answer rows are stored separately:

| File | Queries | Positive answers |
|---|---:|---:|
| `data/papyrus_3prop_2pos/all.csv` | 13,186 | 26,372 |
| `data/papyrus_3prop_2pos/train.csv` | 9,798 | 19,596 |
| `data/papyrus_3prop_2pos/val.csv` | 1,802 | 3,604 |
| `data/papyrus_3prop_2pos/test.csv` | 1,586 | 3,172 |

The Papyrus split reports zero molecule overlap between train/val/test according to `data/papyrus_3prop_2pos/dataset_stats.json`.

Local processing completed:

- PubChem BioAssay archive index exported to CSV.
- PubChem BioAssay full crawl completed: 1,933/1,933 CSV/Data shards downloaded.
- PubChem full-crawl supported observations merged/deduplicated before building `data/pubchem_3prop_2pos`.
- Papyrus++ v05.5 high-quality rows normalized and built into `data/papyrus_3prop_2pos`.
- PubChem compound SDF archive index exported to CSV.
- Tox21 public assay list downloaded and exported to CSV.
- Tox21 public page links exported to CSV.
- ToxCast metadata workbooks converted to CSV.

Not yet included as final rows:

- PubChem BioAssay observations are normalized separately and have a PubChem-only 3-property/2-positive dataset, but they are not merged into `chembl_3prop_2pos` because endpoints are assay-local and still need target/assay enrichment before being treated as the same quality layer as ChEMBL.
- Tox21 assay result tables were not downloaded and normalized.
- ToxCast metadata does not include the complete molecule-level activity matrix with structures.
- DrugBank/eTOX require restricted access.

## Remaining Work To Grow Dataset

The active ChEMBL-side dataset is now in the requested format. PubChem, Papyrus, and BindingDB have separate source-specific outputs; further growth requires target/assay enrichment for these outputs or importing additional molecule-level labels from Tox21 result tables, ToxCast activity matrices, DrugBank/eTOX, or other experimental sources.

Meeting follow-up tracking:

- `MEETING_ACTION_ITEMS.md` lists the requirements from the latest meeting and the acceptance checks for the next report.
- `scripts/generate_meeting_followup_report.py` generates `reports/meeting_followup_report.md` and `reports/meeting_followup_report.json`.
- The first baseline report is available for the BindingDB target-conditioned triplet task. A train-only retrieval/MMP baseline is still needed for the ChEMBL/PubChem/Papyrus 3-property editing tasks.

## BindingDB Target-Conditioned Direction

Status: full BindingDB 202605 TSV snapshot processed into a separate direction-2 dataset.

BindingDB is not merged into the ChEMBL/PubChem 3-property molecular-property editing files. It is stored as a target-conditioned binding/activity ranking dataset:

```text
input ligand SMILES + target protein context + measurement instruction
-> stronger positive ligand + weaker hard-negative ligand
```

Main output folder:

```text
data/bindingdb_target_conditioned/
```

Full-snapshot processing summary:

| Stage | Count |
|---|---:|
| Raw BindingDB rows read | 3,176,528 |
| Parsed Ki/Kd/IC50/EC50 observations | 3,179,049 |
| Exact observations used for ranking | 2,543,532 |
| Median-aggregated ligand-target-measurement rows | 1,853,059 |
| Rank-ready scaffold-filtered rows | 1,807,525 |
| Final ranking pairs | 435,980 |
| Final triplets | 435,980 |

Split summary:

| File | Rows |
|---|---:|
| `train.csv` | 279,525 |
| `val.csv` | 34,019 |
| `test_seen_target.csv` | 34,574 |
| `test_unseen_target.csv` | 46,377 |

Leakage controls in `dataset_stats.json` report zero ligand overlap across all splits and zero target overlap between the seen-target splits and `test_unseen_target`.

Measurement policy:

- `Ki` and `Kd` are grouped as `binding_affinity`.
- `IC50` and `EC50` are grouped as `activity_potency`.
- All nM measurements are converted to p-scale by `p = 9 - log10(value_nM)`.
- Duplicates for the same ligand + target + measurement type are aggregated by median p-scale.
- Censored values such as `<`, `>`, `<=`, `>=` are retained in observations but excluded from pair/triplet ranking.
- BindingDB full TSV does not provide a reliable agonist/antagonist label column, so `modulation_numeric = 0.5` is stored as unknown, not as a measured neutral class.


