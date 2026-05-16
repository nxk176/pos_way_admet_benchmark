#!/usr/bin/env bash
set -euo pipefail

source .venv_lora/bin/activate

# ADMET edit dataset: ChEMBL + PubChem + Papyrus, leakage-filtered after merge.
python scripts/prepare_llm_sft_data.py \
  --preset edit-merged \
  --out-dir data/llm_sft

# BindingDB full SFT dataset: target-conditioned ligand optimization.
python scripts/prepare_llm_sft_data.py \
  --preset bindingdb \
  --out-dir data/llm_sft_bindingdb_full

cat data/llm_sft/admet_edit_3prop_merged/dataset_stats.json
cat data/llm_sft_bindingdb_full/bindingdb_target_conditioned/dataset_stats.json

python server_training/audit_sft_datasets.py \
  --out reports/server_training/dataset_audit.json
