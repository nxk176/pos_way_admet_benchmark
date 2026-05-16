# Data Transfer Manifest

GitHub does not include local data/model artifacts. Transfer these separately to the server.

## Required For Training

The server training workflow needs these folders:

```text
data/chembl_3prop_2pos/
data/pubchem_3prop_2pos/
data/papyrus_3prop_2pos/
data/bindingdb_target_conditioned/
```

These are enough for:

```bash
python scripts/prepare_llm_sft_data.py --preset edit-merged --out-dir data/llm_sft
python scripts/prepare_llm_sft_data.py --preset bindingdb --out-dir data/llm_sft_bindingdb_full
```

## Optional / Reproducibility Data

These are useful for full audit/rebuild, but not required for LoRA training if the final processed datasets above are already present:

```text
raw/
data/bindingdb/
data/pubchem_bioassay/
data/normalized/
data/normalized_csv_expanded/
data/multiproperty_large_rebuild/
data/broad/
data/experimental/
data/properties/
```

## Not Required On Server

GGUF files are not needed for LoRA/QLoRA training:

```text
models/*.gguf
```

The server downloads HF safetensors with:

```bash
python server_training/download_hf_models.py --model all
```
