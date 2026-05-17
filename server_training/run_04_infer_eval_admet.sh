#!/usr/bin/env bash
set -euo pipefail

if [ -f .venv_lora/bin/activate ]; then
  source .venv_lora/bin/activate
else
  echo "No .venv_lora found; using the currently active Python environment."
fi

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export TOKENIZERS_PARALLELISM=false

MAX_INFER_ROWS="${MAX_INFER_ROWS:-0}"
MAX_INFER_ROWS_PER_SOURCE="${MAX_INFER_ROWS_PER_SOURCE:-0}"
MAX_SEQ_LENGTH="${MAX_SEQ_LENGTH:-2048}"
MAX_NEW_TOKENS="${MAX_NEW_TOKENS:-192}"

QWEN14="models/hf/Qwen2.5-14B-Instruct"
DEEPSEEK14="models/hf/DeepSeek-R1-Distill-Qwen-14B"
ADMET_TEST="data/llm_sft/admet_edit_3prop_merged/test.jsonl"

mkdir -p reports/server_training/predictions reports/server_training/metrics

infer_eval () {
  local base_model="$1"
  local adapter="$2"
  local tag="$3"

  local pred="reports/server_training/predictions/${tag}.jsonl"
  local metrics="reports/server_training/metrics/${tag}.metrics.json"

  python server_training/infer_lora.py \
    --base-model "$base_model" \
    --adapter-dir "$adapter" \
    --input-jsonl "$ADMET_TEST" \
    --output-jsonl "$pred" \
    --max-rows "$MAX_INFER_ROWS" \
    --max-rows-per-source "$MAX_INFER_ROWS_PER_SOURCE" \
    --max-seq-length "$MAX_SEQ_LENGTH" \
    --max-new-tokens "$MAX_NEW_TOKENS" \
    --load-in-4bit \
    --bf16 \
    --temperature 0

  python server_training/evaluate_lora_predictions.py \
    --predictions "$pred" \
    --reference-jsonl "$ADMET_TEST" \
    --out "$metrics"
}

infer_eval "$QWEN14" "models/lora/qwen2.5-14b-admet-edit" \
  "qwen2.5-14b-admet-edit_test"

infer_eval "$DEEPSEEK14" "models/lora/deepseek-r1-distill-qwen-14b-admet-edit" \
  "deepseek14b-admet-edit_test"
