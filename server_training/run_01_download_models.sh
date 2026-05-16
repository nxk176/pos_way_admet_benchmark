#!/usr/bin/env bash
set -euo pipefail

source .venv_lora/bin/activate

# Optional: export HF_TOKEN=... if your server needs authenticated HF downloads.
python server_training/download_hf_models.py --model all

du -sh models/hf/* || true
