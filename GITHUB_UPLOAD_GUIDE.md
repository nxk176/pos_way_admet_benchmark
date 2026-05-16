# GitHub Upload Guide

This repository contains large local datasets and model files. GitHub should track code, documentation, schemas, and reproducible scripts only.

Do not push these folders to GitHub:

```text
data/
raw/
models/
reports/
```

They are ignored by `.gitignore`.

## What To Upload To GitHub

Track:

```text
configs/
schemas/
scripts/
server_training/
README.md
*.md
```

Do not track:

```text
*.gguf
*.safetensors
*.zip
*.csv
*.jsonl
*.db
```

## One-Time GitHub Push

From this folder:

```powershell
cd C:\Users\ADMIN\Desktop\rebuild\pos_way_admet_benchmark
git init
git branch -M main
git add .
git status
git commit -m "Initial POS-WAY benchmark workflow"
git remote add origin https://github.com/YOUR_USERNAME/YOUR_REPO_NAME.git
git push -u origin main
```

Before `git commit`, inspect `git status`. If `data/`, `raw/`, `models/`, or `reports/` appear as staged files, stop and fix `.gitignore`.

## Clone On Server

```bash
git clone https://github.com/YOUR_USERNAME/YOUR_REPO_NAME.git
cd YOUR_REPO_NAME
```

Then transfer the local data folder separately, for example:

```bash
scp -r C:/Users/ADMIN/Desktop/rebuild/pos_way_admet_benchmark/data user@server:/path/to/YOUR_REPO_NAME/
```

Or use `rsync` from Linux/WSL:

```bash
rsync -av --progress pos_way_admet_benchmark/data/ user@server:/path/to/YOUR_REPO_NAME/data/
```

Model safetensors do not need to be uploaded. On the server, run:

```bash
python server_training/download_hf_models.py --model all
```

## Server Run Order

```bash
cd YOUR_REPO_NAME
python -m pip install -r server_training/requirements.txt
python server_training/download_hf_models.py --model all
python scripts/prepare_llm_sft_data.py --preset edit-merged --out-dir data/llm_sft
python scripts/prepare_llm_sft_data.py --preset bindingdb --out-dir data/llm_sft_bindingdb_full
python server_training/audit_sft_datasets.py --out reports/server_training/dataset_audit.json
```

Then train:

```bash
MAX_TRAIN_SAMPLES=2000 MAX_VAL_SAMPLES=200 bash server_training/run_03_train_all.sh
MAX_INFER_ROWS=100 bash server_training/run_04_infer_eval_all.sh
```

If the pilot run is clean:

```bash
bash server_training/run_03_train_all.sh
bash server_training/run_04_infer_eval_all.sh
```
