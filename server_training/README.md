# Server Training Workflow

This folder is the server-ready workflow for LoRA/QLoRA experiments.

It trains two datasets:

1. `admet_edit_3prop_merged`: ChEMBL + PubChem + Papyrus molecular editing.
2. `bindingdb_target_conditioned`: BindingDB target-conditioned ligand optimization.

It fine-tunes two HF safetensors models:

1. `Qwen/Qwen2.5-14B-Instruct`
2. `deepseek-ai/DeepSeek-R1-Distill-Qwen-14B`

GGUF files are not used for training. GGUF is only for local inference with llama.cpp / LM Studio / Ollama.

## Expected Server

Recommended:

```text
NVIDIA A30 24GB VRAM
Linux
CUDA-compatible PyTorch
Python 3.10+
```

## Model Sources

The scripts download:

- `Qwen/Qwen2.5-14B-Instruct`
- `deepseek-ai/DeepSeek-R1-Distill-Qwen-14B`

These are the HF safetensors checkpoints needed for LoRA/QLoRA.

## Quick Start

From repo root:

```bash
cd pos_way_admet_benchmark
bash server_training/run_00_setup_env.sh
bash server_training/run_01_download_models.sh
bash server_training/run_02_prepare_datasets.sh
bash server_training/run_03_train_all.sh
bash server_training/run_04_infer_eval_all.sh
```

If Hugging Face requires authentication:

```bash
export HF_TOKEN=your_token_here
huggingface-cli login
```

## Outputs

Model adapters:

```text
models/lora/qwen2.5-14b-admet-edit/
models/lora/deepseek-r1-distill-qwen-14b-admet-edit/
models/lora/qwen2.5-14b-bindingdb/
models/lora/deepseek-r1-distill-qwen-14b-bindingdb/
```

Predictions:

```text
reports/server_training/predictions/
```

Metrics:

```text
reports/server_training/metrics/
```

## Dataset Preparation

`run_02_prepare_datasets.sh` writes:

```text
data/llm_sft/admet_edit_3prop_merged/
data/llm_sft_bindingdb_full/bindingdb_target_conditioned/
```

ADMET merge policy:

```text
ChEMBL + PubChem + Papyrus are merged because they share:
SMILES + instruction -> exactly 2 edited SMILES
```

BindingDB is kept separate because it has target protein context:

```text
ligand + protein target -> stronger ligand + hard negative
```

`run_02_prepare_datasets.sh` also writes:

```text
reports/server_training/dataset_audit.json
```

The merged ADMET training set is intentionally allowed to use all valid source rows, but the test set is source-imbalanced because PubChem is much larger than ChEMBL/Papyrus. For paper reporting, use both:

- micro-average metrics over the merged test set
- per-source metrics and source-macro metrics

## Pilot Runs

Before a full run, use a short pilot:

```bash
MAX_TRAIN_SAMPLES=2000 MAX_VAL_SAMPLES=200 MAX_INFER_ROWS=100 \
bash server_training/run_03_train_all.sh
```

Then:

```bash
MAX_INFER_ROWS=100 bash server_training/run_04_infer_eval_all.sh
```

For the full paper run, leave `MAX_TRAIN_SAMPLES=0`, which means no training subsample.

## Single Experiment Example

Qwen 14B on ADMET:

```bash
source .venv_lora/bin/activate

python server_training/train_qlora.py \
  --base-model models/hf/Qwen2.5-14B-Instruct \
  --train-jsonl data/llm_sft/admet_edit_3prop_merged/train.jsonl \
  --val-jsonl data/llm_sft/admet_edit_3prop_merged/val.jsonl \
  --output-dir models/lora/qwen2.5-14b-admet-edit \
  --load-in-4bit \
  --gradient-checkpointing \
  --bf16 \
  --max-seq-length 2048 \
  --num-train-epochs 1 \
  --per-device-train-batch-size 1 \
  --gradient-accumulation-steps 16 \
  --learning-rate 2e-4 \
  --lora-r 16 \
  --lora-alpha 32
```

Inference:

```bash
python server_training/infer_lora.py \
  --base-model models/hf/Qwen2.5-14B-Instruct \
  --adapter-dir models/lora/qwen2.5-14b-admet-edit \
  --input-jsonl data/llm_sft/admet_edit_3prop_merged/test.jsonl \
  --output-jsonl reports/server_training/predictions/qwen2.5-14b-admet-edit_test.jsonl \
  --max-rows 1000 \
  --load-in-4bit \
  --bf16 \
  --temperature 0
```

Evaluation:

```bash
python server_training/evaluate_lora_predictions.py \
  --predictions reports/server_training/predictions/qwen2.5-14b-admet-edit_test.jsonl \
  --out reports/server_training/metrics/qwen2.5-14b-admet-edit_test.metrics.json
```

## Notes For Reporting

Report these as separate baselines:

- Qwen2.5-14B-Instruct LoRA on ADMET edit.
- DeepSeek-R1-Distill-Qwen-14B LoRA on ADMET edit.
- Qwen2.5-14B-Instruct LoRA on BindingDB.
- DeepSeek-R1-Distill-Qwen-14B LoRA on BindingDB.

Do not merge ADMET and BindingDB into one fine-tuning dataset for the first serious baseline table.
