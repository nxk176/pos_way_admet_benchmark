#!/usr/bin/env bash
set -euo pipefail

python -m venv .venv_lora
source .venv_lora/bin/activate
python -m pip install --upgrade pip
python -m pip install -r server_training/requirements.txt

python - <<'PY'
import torch
print("torch:", torch.__version__)
print("cuda_available:", torch.cuda.is_available())
if torch.cuda.is_available():
    print("gpu:", torch.cuda.get_device_name(0))
    print("capability:", torch.cuda.get_device_capability(0))
PY
