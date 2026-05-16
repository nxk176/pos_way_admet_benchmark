# POS-WAY 3-Property 2-Positive Multi-Answer Dataset

Current canonical schema:

```text
one input molecule
+ one English instruction with exactly 3 property objective groups
-> exactly 2 positive answer molecules
```

The negative candidate requirement was removed because enforcing `2 positive + 1 negative` dropped the usable query count below the requested scale.

## ChEMBL 3-Property 2-Positive Dataset

| File | Queries | Positive answers | Unique inputs | Primary endpoints | Preserved secondary endpoints |
|---|---:|---:|---:|---:|---:|
| `data/chembl_3prop_2pos/all.csv` | 6,070 | 12,140 | 3,994 | 237 | 261 |
| `data/chembl_3prop_2pos/train.csv` | 5,858 | 11,716 | 3,803 | 235 | 258 |
| `data/chembl_3prop_2pos/val.csv` | 104 | 208 | 99 | 23 | 22 |
| `data/chembl_3prop_2pos/test.csv` | 108 | 216 | 92 | 29 | 28 |
| `data/chembl_3prop_2pos/dataset_stats.json` | - | - | - | - | - |

Each query has exactly two positive answers, so `all.csv` has 12,140 answer slots.
The loose pretrain-style ChEMBL rows are folded into `train.csv`; `val.csv` and `test.csv` remain strict held-out splits with zero molecule overlap against train.

## Row Format

Main columns:

| Column | Meaning |
|---|---|
| `query_id` | Unique query ID. |
| `split` | `train`, `val`, or `test`; `all.csv` is the union of those splits. |
| `instruction` | English multi-property instruction. |
| `input_smiles_canon` | Input molecule. |
| `primary_endpoint` | Experimental property to improve. |
| `primary_objective_json` | Primary endpoint objective and measured delta evidence. |
| `preserved_property_json` | One experimental secondary property to preserve. |
| `local_constraints_json` | Local constraints: MW, LogP, QED, and SA. |
| `num_property_objectives` | Always `3`. |
| `num_positive_answers` | Always `2`. |
| `positive_answer_smiles_json` | JSON list with exactly two accepted output SMILES. |
| `positive_answers_json` | Full evidence for the two answer molecules. |

Instruction template:

```text
Increase/Decrease {primary_endpoint} while preserving {secondary_endpoint} within {tolerance}; also keep MW, LogP, QED, and synthetic accessibility within local edit constraints.
```

## Remaining Sources Folder

Additional non-ZINC sources that can be locally processed are stored separately:

```text
data/remaining_sources_3prop_2pos/
```

That folder contains PubChem index tables, Tox21 assay list/page links, and converted ToxCast metadata workbooks.

PubChem BioAssay also has a separate normalized observation output:

| File | Rows | Meaning |
|---|---:|---|
| `data/pubchem_bioassay/pubchem_bioassay_observations.csv` | 127,360 | Molecule-level PubChem BioAssay observations from downloaded CSV/Data shards. |
| `data/pubchem_bioassay/pubchem_property_summary.csv` | 27,831 | Assay-local PubChem endpoint summary table. |
| `data/pubchem_bioassay/pubchem_final_supported_observations_merged.csv` | 2,879,222 | Full-crawl exact p-scale PubChem observations used for the PubChem source-specific dataset. |

These PubChem rows are not merged into `chembl_3prop_2pos`; they are kept as a separate source layer because assay-local endpoints still need target/assay enrichment before being treated as the same quality layer as ChEMBL/BindingDB.

## PubChem 3-Property 2-Positive Dataset

PubChem has a separate source-specific view for inspection:

| File | Queries | Positive answers | Unique inputs | Primary endpoints | Preserved secondary endpoints |
|---|---:|---:|---:|---:|---:|
| `data/pubchem_3prop_2pos/all.csv` | 42,754 | 85,508 | 19,256 | 2,544 | 3,266 |
| `data/pubchem_3prop_2pos/train.csv` | 34,204 | 68,408 | 15,547 | 2,259 | 2,790 |
| `data/pubchem_3prop_2pos/val.csv` | 4,275 | 8,550 | 1,830 | 470 | 536 |
| `data/pubchem_3prop_2pos/test.csv` | 4,275 | 8,550 | 1,879 | 460 | 509 |
| `data/pubchem_3prop_2pos/dataset_stats.json` | - | - | - | - | - |

This view uses the same `1 input + 3-property instruction -> exactly 2 positive answers` shape, but it is PubChem-only and remains separate from the active ChEMBL-side dataset.

## Papyrus 3-Property 2-Positive Dataset

Papyrus has a separate source-specific view for inspection:

| File | Queries | Positive answers | Unique inputs | Primary endpoints | Preserved secondary endpoints |
|---|---:|---:|---:|---:|---:|
| `data/papyrus_3prop_2pos/all.csv` | 13,186 | 26,372 | 7,540 | 887 | 13,186 |
| `data/papyrus_3prop_2pos/train.csv` | 9,798 | 19,596 | - | - | - |
| `data/papyrus_3prop_2pos/val.csv` | 1,802 | 3,604 | - | - | - |
| `data/papyrus_3prop_2pos/test.csv` | 1,586 | 3,172 | - | - | - |
| `data/papyrus_3prop_2pos/dataset_stats.json` | - | - | - | - | - |

This view uses the same `1 input + 3-property instruction -> exactly 2 positive answers` shape, but it is Papyrus-only and remains separate from the active ChEMBL-side dataset.

## Meeting Follow-up Report

The latest meeting follow-up can be regenerated with:

```powershell
python .\scripts\generate_meeting_followup_report.py
```

It writes `reports/meeting_followup_report.md` and `reports/meeting_followup_report.json`, including dataset inventory, shape audits, and first BindingDB heuristic baselines.

## Traditional Baseline Training

Traditional candidate-ranking baselines are documented in `BASELINE_GUIDE.md` and implemented in:

```text
scripts/train_traditional_baselines.py
```

The script supports ChEMBL/PubChem/Papyrus 3-property edit ranking and BindingDB target-conditioned pairwise ranking. It is intended for supervised/classical baselines using molecule graph fingerprints plus instruction/metadata features, not de novo generation.

## LLM Direct-Edit Baselines

LLM direct-edit baselines are documented in `LLM_BASELINE_GUIDE.md` and implemented in:

```text
scripts/run_llm_edit_baseline.py
```

This script evaluates the task shape expected for GPT/Qwen-style baselines:

```text
input SMILES + natural-language instruction -> exactly 2 edited SMILES
```

It supports prompt export, OpenAI-compatible API calls, local HuggingFace/Qwen inference, and evaluation with validity, exact-hit, Tanimoto-to-gold, and local-constraint metrics.

Local GGUF model choices, current downloaded files, and Hugging Face links are tracked in `LOCAL_LLM_MODELS.md`.

Step-by-step local GGUF runs, dataset merge policy, SFT JSONL preparation, and LoRA commands are documented in `LLM_RUN_AND_FINETUNE_GUIDE.md`.

Server-ready A30/24GB LoRA workflow files are in `server_training/`, including model download, dataset preparation, four training runs, inference, evaluation, and `requirements.txt`.

## Source Tables Used

| Source layer | Path | Rows |
|---|---|---:|
| Expanded normalized observations | `data/normalized_csv_expanded/property_observations.csv` | 697,404 |
| Expanded property summary | `data/normalized_csv_expanded/property_summary.csv` | 8,232 endpoints |
| ChEMBL broad pChEMBL | `data/broad/chembl_broad_pchembl_observations.csv` | 500,000 |
| BindingDB curated observations | `data/bindingdb/bindingdb_curated_observations.csv` | 95,417 |
| ChEMBL ADMET-like measurements | `data/experimental/chembl_admet_measurements.jsonl` | 78,692 |
| ChEMBL DRD2/GSK3B/JNK3 properties | `data/properties/chembl_target_activity_properties.jsonl` | 24,972 |
| PubChem BioAssay normalized observations | `data/pubchem_bioassay/pubchem_bioassay_observations.csv` | 127,360 |
| PubChem full-crawl supported observations | `data/pubchem_bioassay/pubchem_final_supported_observations_merged.csv` | 2,879,222 |
| PubChem source-specific 3-property/2-positive rows | `data/pubchem_3prop_2pos/all.csv` | 42,754 queries |
| Papyrus source-specific 3-property/2-positive rows | `data/papyrus_3prop_2pos/all.csv` | 13,186 queries |
| BindingDB target-conditioned triplets | `data/bindingdb_target_conditioned/all.csv` | 435,980 triplets |

## Schema

Use `schemas/chembl_3prop_2pos.schema.json` for one CSV row.



