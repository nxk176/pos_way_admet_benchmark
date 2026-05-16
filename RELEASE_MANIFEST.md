# Release Manifest

## Active ChEMBL 3-Property 2-Positive Dataset

Schema:

```text
1 input molecule + 3-property instruction -> exactly 2 positive answers
```

| File | Queries | Positive answers |
|---|---:|---:|
| `data/chembl_3prop_2pos/all.csv` | 6,070 | 12,140 |
| `data/chembl_3prop_2pos/train.csv` | 5,858 | 11,716 |
| `data/chembl_3prop_2pos/val.csv` | 104 | 208 |
| `data/chembl_3prop_2pos/test.csv` | 108 | 216 |
| `data/chembl_3prop_2pos/dataset_stats.json` | - | - |

Schema file: `schemas/chembl_3prop_2pos.schema.json`.

## Remaining Sources Folder

`data/remaining_sources_3prop_2pos/` contains processed non-ZINC source artifacts that are not yet final rows:

- PubChem BioAssay CSV archive index.
- PubChem compound SDF archive index.
- Tox21 public assay list and page links.
- ToxCast metadata workbooks converted to CSV.
- `source_processing_status.csv` explaining why each source is or is not usable as final rows right now.

## PubChem Separate 3-Property 2-Positive Dataset

PubChem has been processed into its own source-specific dataset for review, without merging it into the active ChEMBL-side set:

| File | Queries | Positive answers |
|---|---:|---:|
| `data/pubchem_3prop_2pos/all.csv` | 42,754 | 85,508 |
| `data/pubchem_3prop_2pos/train.csv` | 34,204 | 68,408 |
| `data/pubchem_3prop_2pos/val.csv` | 4,275 | 8,550 |
| `data/pubchem_3prop_2pos/test.csv` | 4,275 | 8,550 |
| `data/pubchem_3prop_2pos/dataset_stats.json` | - | - |

Source observation layer:

- `data/pubchem_bioassay/pubchem_bioassay_observations.csv`
- `data/pubchem_bioassay/pubchem_property_summary.csv`
- `data/pubchem_bioassay/pubchem_final_supported_observations_merged.csv`

PubChem crawl coverage: 1,933/1,933 BioAssay CSV/Data shards downloaded. The full-crawl supported observation file has 2,879,222 deduplicated exact p-scale rows.

## Papyrus Separate 3-Property 2-Positive Dataset

Papyrus has been processed into its own source-specific dataset for review, without merging it into the active ChEMBL-side set:

| File | Queries | Positive answers |
|---|---:|---:|
| `data/papyrus_3prop_2pos/all.csv` | 13,186 | 26,372 |
| `data/papyrus_3prop_2pos/train.csv` | 9,798 | 19,596 |
| `data/papyrus_3prop_2pos/val.csv` | 1,802 | 3,604 |
| `data/papyrus_3prop_2pos/test.csv` | 1,586 | 3,172 |
| `data/papyrus_3prop_2pos/dataset_stats.json` | - | - |

Source snapshot: Papyrus++ high-quality v05.5. The split leakage check reports zero molecule overlap.

## BindingDB Target-Conditioned Dataset

BindingDB is released as a separate target-conditioned ranking/triplet dataset:

| File | Rows |
|---|---:|
| `data/bindingdb_target_conditioned/all.csv` | 435,980 |
| `data/bindingdb_target_conditioned/train.csv` | 279,525 |
| `data/bindingdb_target_conditioned/val.csv` | 34,019 |
| `data/bindingdb_target_conditioned/test_seen_target.csv` | 34,574 |
| `data/bindingdb_target_conditioned/test_unseen_target.csv` | 46,377 |
| `data/bindingdb_target_conditioned/dataset_stats.json` | - |

## Source/Audit Layers

- `data/normalized_csv_expanded/property_observations.csv`
- `data/normalized_csv_expanded/property_summary.csv`
- `data/broad/chembl_broad_pchembl_observations.csv`
- `data/bindingdb/bindingdb_curated_observations.csv`
- `data/experimental/chembl_admet_measurements.jsonl`
- `data/properties/chembl_target_activity_properties.jsonl`
- `raw/public/`
- `configs/`
- `scripts/`

## Documentation

- `README.md`
- `SOURCE_DATA_FILES.md`
- `REQUIREMENT_STATUS.md`
- `PROPERTY_SUMMARY_TABLE_EXPANDED.html`
- `MEETING_ACTION_ITEMS.md`
- `BASELINE_GUIDE.md`
- `LLM_BASELINE_GUIDE.md`
- `LLM_RUN_AND_FINETUNE_GUIDE.md`
- `LOCAL_LLM_MODELS.md`
- `server_training/`
- `reports/meeting_followup_report.md`



