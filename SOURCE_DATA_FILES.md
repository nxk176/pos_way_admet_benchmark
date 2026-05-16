# Source-Specific Data Files

The active ChEMBL-side dataset is:

```text
1 input molecule + 3-property instruction -> exactly 2 positive answers
```

## ChEMBL 3-Property 2-Positive Dataset

| File | Queries | Positive answers | Function |
|---|---:|---:|---|
| `data/chembl_3prop_2pos/all.csv` | 6,070 | 12,140 | ChEMBL-only union of train/val/test. |
| `data/chembl_3prop_2pos/train.csv` | 5,858 | 11,716 | Train split: strict train plus leakage-filtered loose pretrain-style rows. |
| `data/chembl_3prop_2pos/val.csv` | 104 | 208 | Strict validation split. |
| `data/chembl_3prop_2pos/test.csv` | 108 | 216 | Strict held-out split. |
| `data/chembl_3prop_2pos/dataset_stats.json` | - | - | Counts and build policy. |

Every row has:

- exactly 3 objective groups: primary experimental property, one preserved experimental property, local MW/LogP/QED/SA constraints;
- exactly 2 positive answer molecules;
- no negative candidate output.

## Remaining Sources

| Folder/File | Content |
|---|---|
| `data/remaining_sources_3prop_2pos/source_processing_status.csv` | Per-source status for PubChem, Tox21, ToxCast, DrugBank, eTOX. |
| `data/remaining_sources_3prop_2pos/pubchem_bioassay_csv_index.csv` | PubChem BioAssay archive index. |
| `data/remaining_sources_3prop_2pos/pubchem_compound_sdf_index.csv` | PubChem compound SDF archive index. |
| `data/remaining_sources_3prop_2pos/tox21_assays.csv` | Tox21 public assay list downloaded from the public app API. |
| `data/remaining_sources_3prop_2pos/tox21_public_page_links.csv` | Links extracted from the Tox21 public page. |
| `data/remaining_sources_3prop_2pos/toxcast_metadata_csv/` | ToxCast metadata workbooks converted to CSV. |

These remaining-source files are stored separately because they are not yet molecule-level edit labels.

## PubChem Separate 3-Property 2-Positive Dataset

This is not merged with the active ChEMBL-side dataset, but it has the same row shape for separate review:

| File | Queries | Positive answers | Function |
|---|---:|---:|---|
| `data/pubchem_3prop_2pos/all.csv` | 42,754 | 85,508 | PubChem-only 3-property/2-positive view. |
| `data/pubchem_3prop_2pos/train.csv` | 34,204 | 68,408 | PubChem train split. |
| `data/pubchem_3prop_2pos/val.csv` | 4,275 | 8,550 | PubChem validation split. |
| `data/pubchem_3prop_2pos/test.csv` | 4,275 | 8,550 | PubChem held-out split. |
| `data/pubchem_3prop_2pos/dataset_stats.json` | - | - | Counts, policy, and leakage check. |

PubChem selection rule: exact PubChem p-scale primary measurements, same assay-local endpoint, same BRICS/MMP core, one preserved PubChem secondary endpoint, and MW/LogP/QED/SA local constraints.
The PubChem source was full-crawled from 1,933/1,933 BioAssay CSV/Data shards and merged/deduplicated before building this view.

## Papyrus Separate 3-Property 2-Positive Dataset

This is not merged with the active ChEMBL-side dataset, but it has the same row shape for separate review:

| File | Queries | Positive answers | Function |
|---|---:|---:|---|
| `data/papyrus_3prop_2pos/all.csv` | 13,186 | 26,372 | Papyrus-only 3-property/2-positive view. |
| `data/papyrus_3prop_2pos/train.csv` | 9,798 | 19,596 | Papyrus train split. |
| `data/papyrus_3prop_2pos/val.csv` | 1,802 | 3,604 | Papyrus validation split. |
| `data/papyrus_3prop_2pos/test.csv` | 1,586 | 3,172 | Papyrus held-out split. |
| `data/papyrus_3prop_2pos/dataset_stats.json` | - | - | Counts, policy, and leakage check. |

Papyrus selection rule: exact quantitative Papyrus++ high-quality bioactivity records, same target/assay-derived endpoint, shared Murcko scaffold, one preserved secondary endpoint, and MW/LogP/QED/SA local constraints.

## BindingDB Target-Conditioned Dataset

BindingDB is a separate target-conditioned ranking task, not a 3-property editing table:

```text
input ligand + target protein context + measurement instruction -> stronger positive ligand + weaker hard-negative ligand
```

| File | Rows | Function |
|---|---:|---|
| `data/bindingdb_target_conditioned/all.csv` | 435,980 | Full triplet table. |
| `data/bindingdb_target_conditioned/train.csv` | 279,525 | Train split. |
| `data/bindingdb_target_conditioned/val.csv` | 34,019 | Validation split. |
| `data/bindingdb_target_conditioned/test_seen_target.csv` | 34,574 | Held-out ligands on seen targets. |
| `data/bindingdb_target_conditioned/test_unseen_target.csv` | 46,377 | Held-out ligands on unseen targets. |
| `data/bindingdb_target_conditioned/dataset_stats.json` | - | Build stats, leakage checks, and measurement policy. |

## Source/Audit Files

| File | Content |
|---|---|
| `data/normalized_csv_expanded/property_observations.csv` | 697,404 normalized observations. |
| `data/normalized_csv_expanded/property_summary.csv` | Summary of 8,232 endpoints. |
| `PROPERTY_SUMMARY_TABLE_EXPANDED.html` | Browser-friendly endpoint summary. |
| `data/broad/chembl_broad_pchembl_observations.csv` | ChEMBL broad pChEMBL source layer. |
| `data/bindingdb/bindingdb_curated_observations.csv` | BindingDB curated source layer. |
| `data/experimental/chembl_admet_measurements.jsonl` | ChEMBL ADMET-like source/audit layer. |
| `data/properties/chembl_target_activity_properties.jsonl` | ChEMBL DRD2/GSK3B/JNK3 source/audit layer. |
| `data/pubchem_bioassay/pubchem_final_supported_observations_merged.csv` | 2,879,222 full-crawl PubChem exact p-scale observations used for the PubChem source-specific view. |
| `data/papyrus_3prop_2pos/dataset_stats.json` | Papyrus source-specific build stats. |
| `reports/meeting_followup_report.md` | Meeting follow-up inventory, audit, and first BindingDB baseline report. |



