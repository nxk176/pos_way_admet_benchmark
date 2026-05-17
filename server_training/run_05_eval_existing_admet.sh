#!/usr/bin/env bash
set -euo pipefail

if [ -f .venv_lora/bin/activate ]; then
  source .venv_lora/bin/activate
else
  echo "No .venv_lora found; using the currently active Python environment."
fi

PRED_DIR="${PRED_DIR:-reports/server_training/predictions}"
METRIC_DIR="${METRIC_DIR:-reports/server_training/metrics}"
ADMET_TEST="${ADMET_TEST:-data/llm_sft/admet_edit_3prop_merged/test.jsonl}"

mkdir -p "$METRIC_DIR"

evaluate_one () {
  local tag="$1"
  local pred="$PRED_DIR/${tag}.jsonl"
  local metrics="$METRIC_DIR/${tag}.metrics.json"

  if [ ! -f "$pred" ]; then
    echo "ERROR: missing prediction file: $pred" >&2
    return 1
  fi

  python server_training/evaluate_lora_predictions.py \
    --predictions "$pred" \
    --reference-jsonl "$ADMET_TEST" \
    --out "$metrics"
}

evaluate_one "qwen2.5-14b-admet-edit_test"
evaluate_one "deepseek14b-admet-edit_test"
