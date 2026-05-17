#!/usr/bin/env bash
set -euo pipefail

if [ -f .venv_lora/bin/activate ]; then
  source .venv_lora/bin/activate
else
  echo "No .venv_lora found; using the currently active Python environment."
fi

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export TOKENIZERS_PARALLELISM=false
export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"

MAX_TRAIN_SAMPLES="${MAX_TRAIN_SAMPLES:-0}"
MAX_VAL_SAMPLES="${MAX_VAL_SAMPLES:-2000}"
EPOCHS="${EPOCHS:-1}"
MAX_SEQ_LENGTH="${MAX_SEQ_LENGTH:-2048}"
GRAD_ACCUM="${GRAD_ACCUM:-16}"
LR="${LR:-2e-4}"

QWEN14="models/hf/Qwen2.5-14B-Instruct"
DEEPSEEK14="models/hf/DeepSeek-R1-Distill-Qwen-14B"

ADMET_TRAIN="data/llm_sft/admet_edit_3prop_merged/train.jsonl"
ADMET_VAL="data/llm_sft/admet_edit_3prop_merged/val.jsonl"

train_one () {
  local base_model="$1"
  local output_dir="$2"
  local run_name="$3"

  python server_training/train_qlora.py \
    --base-model "$base_model" \
    --train-jsonl "$ADMET_TRAIN" \
    --val-jsonl "$ADMET_VAL" \
    --output-dir "$output_dir" \
    --run-name "$run_name" \
    --load-in-4bit \
    --gradient-checkpointing \
    --bf16 \
    --max-seq-length "$MAX_SEQ_LENGTH" \
    --num-train-epochs "$EPOCHS" \
    --per-device-train-batch-size 1 \
    --per-device-eval-batch-size 1 \
    --gradient-accumulation-steps "$GRAD_ACCUM" \
    --learning-rate "$LR" \
    --lora-r 16 \
    --lora-alpha 32 \
    --lora-dropout 0.05 \
    --eval-steps 250 \
    --save-steps 250 \
    --save-total-limit 2 \
    --max-train-samples "$MAX_TRAIN_SAMPLES" \
    --max-val-samples "$MAX_VAL_SAMPLES"
}

train_one "$QWEN14" \
  "models/lora/qwen2.5-14b-admet-edit" "qwen2.5-14b-admet-edit"

train_one "$DEEPSEEK14" \
  "models/lora/deepseek-r1-distill-qwen-14b-admet-edit" "deepseek14b-admet-edit"
