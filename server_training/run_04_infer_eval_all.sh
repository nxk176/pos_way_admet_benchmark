#!/usr/bin/env bash
set -euo pipefail

if [ -f .venv_lora/bin/activate ]; then
  source .venv_lora/bin/activate
else
  echo "No .venv_lora found; using the currently active Python environment."
fi

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
export TOKENIZERS_PARALLELISM=false

MAX_INFER_ROWS="${MAX_INFER_ROWS:-1000}"
MAX_SEQ_LENGTH="${MAX_SEQ_LENGTH:-2048}"

QWEN14="models/hf/Qwen2.5-14B-Instruct"
DEEPSEEK14="models/hf/DeepSeek-R1-Distill-Qwen-14B"

ADMET_TEST="data/llm_sft/admet_edit_3prop_merged/test.jsonl"
BIND_SEEN="data/llm_sft_bindingdb_full/bindingdb_target_conditioned/test_seen_target.jsonl"
BIND_UNSEEN="data/llm_sft_bindingdb_full/bindingdb_target_conditioned/test_unseen_target.jsonl"

mkdir -p reports/server_training/predictions reports/server_training/metrics

infer_eval () {
  local base_model="$1"
  local adapter="$2"
  local input_jsonl="$3"
  local tag="$4"

  local pred="reports/server_training/predictions/${tag}.jsonl"
  local metrics="reports/server_training/metrics/${tag}.metrics.json"

  python server_training/infer_lora.py \
    --base-model "$base_model" \
    --adapter-dir "$adapter" \
    --input-jsonl "$input_jsonl" \
    --output-jsonl "$pred" \
    --max-rows "$MAX_INFER_ROWS" \
    --max-seq-length "$MAX_SEQ_LENGTH" \
    --max-new-tokens 192 \
    --load-in-4bit \
    --bf16 \
    --temperature 0

  python server_training/evaluate_lora_predictions.py \
    --predictions "$pred" \
    --out "$metrics"
}

infer_eval "$QWEN14" "models/lora/qwen2.5-14b-admet-edit" \
  "$ADMET_TEST" "qwen2.5-14b-admet-edit_test"

infer_eval "$DEEPSEEK14" "models/lora/deepseek-r1-distill-qwen-14b-admet-edit" \
  "$ADMET_TEST" "deepseek14b-admet-edit_test"

infer_eval "$QWEN14" "models/lora/qwen2.5-14b-bindingdb" \
  "$BIND_SEEN" "qwen2.5-14b-bindingdb_seen"

infer_eval "$QWEN14" "models/lora/qwen2.5-14b-bindingdb" \
  "$BIND_UNSEEN" "qwen2.5-14b-bindingdb_unseen"

infer_eval "$DEEPSEEK14" "models/lora/deepseek-r1-distill-qwen-14b-bindingdb" \
  "$BIND_SEEN" "deepseek14b-bindingdb_seen"

infer_eval "$DEEPSEEK14" "models/lora/deepseek-r1-distill-qwen-14b-bindingdb" \
  "$BIND_UNSEEN" "deepseek14b-bindingdb_unseen"
