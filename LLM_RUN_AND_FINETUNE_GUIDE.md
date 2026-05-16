# LLM Run And Fine-Tune Guide

This guide covers three separate things:

1. Which datasets can be merged.
2. How to run the downloaded GGUF models as zero-shot / few-shot baselines.
3. How to fine-tune with LoRA/QLoRA.

## 1. Dataset Merge Decision

Recommended merge:

| Dataset | Merge decision | Reason |
|---|---|---|
| `data/chembl_3prop_2pos` | Merge | Same task shape: input SMILES + 3-property instruction -> exactly 2 edited SMILES. |
| `data/pubchem_3prop_2pos` | Merge | Same task shape, source retained in metadata. |
| `data/papyrus_3prop_2pos` | Merge | Same task shape, source retained in metadata. |
| `data/bindingdb_target_conditioned` | Keep separate | Different task: ligand + protein target context -> stronger ligand plus hard negative. |

Primary fine-tuning dataset:

```text
data/llm_sft/admet_edit_3prop_merged/
```

Generated split sizes:

| Split | Rows | Notes |
|---|---:|---|
| train | 48,026 | ChEMBL + PubChem + Papyrus, after dropping train rows overlapping held-out molecules. |
| val | 6,181 | Merged validation set. |
| test | 5,910 | Merged test set after dropping 59 rows overlapping validation molecules. |

Leakage check after merge:

```text
train-val overlap: 0
train-test overlap: 0
val-test overlap: 0
```

BindingDB pilot SFT dataset:

```text
data/llm_sft_bindingdb_20k/bindingdb_target_conditioned/
```

| Split | Rows |
|---|---:|
| train | 20,000 |
| val | 2,000 |
| test_seen_target | 2,000 |
| test_unseen_target | 2,000 |

Regenerate the merged ADMET SFT data:

```powershell
cd C:\Users\ADMIN\Desktop\rebuild\pos_way_admet_benchmark
python .\scripts\prepare_llm_sft_data.py --preset edit-merged --out-dir data/llm_sft
```

Regenerate the BindingDB pilot SFT data:

```powershell
python .\scripts\prepare_llm_sft_data.py `
  --preset bindingdb `
  --max-train-per-source 20000 `
  --max-val-per-source 2000 `
  --max-test-per-source 2000 `
  --out-dir data/llm_sft_bindingdb_20k
```

## 2. Current GGUF Model Status

Current files found in `models/`:

| File | Status |
|---|---|
| `DeepSeek-R1-Distill-Qwen-14B-Q3_K_M.gguf` | Ready to host as one GGUF file. |
| `qwen2.5-14b-instruct-q3_k_m-00001-of-00002.gguf` + `00002-of-00002` | Both shards are present. Point llama.cpp to shard `00001-of-00002`. |
| `qwen2.5-7b-instruct-q5_k_m-00001-of-00002.gguf` + `00002-of-00002` | Both shards are present. Point llama.cpp to shard `00001-of-00002`. |

Check shard status:

```powershell
Get-ChildItem .\models -File | Select-Object Name,Length
```

## 3. Host A GGUF Model Locally

`llama-server`, `llama-cli`, and `ollama` were not found in PATH during this check. Install one host first, or use LM Studio with its OpenAI-compatible local server.

Example with `llama-server`:

```powershell
cd C:\Users\ADMIN\Desktop\rebuild\pos_way_admet_benchmark
llama-server -m .\models\DeepSeek-R1-Distill-Qwen-14B-Q3_K_M.gguf -c 4096 --port 8080
```

For a split Qwen GGUF, put both shards in the same folder, then point `-m` to shard `00001-of-00002`.

## 4. Run Zero-Shot / Few-Shot On ADMET Edit Datasets

Open a second terminal after the local server is running:

```powershell
cd C:\Users\ADMIN\Desktop\rebuild\pos_way_admet_benchmark
$env:LOCAL_LLM_API_KEY="dummy"
```

Run a small smoke test first:

```powershell
python .\scripts\run_llm_edit_baseline.py run `
  --provider openai-compatible `
  --base-url http://127.0.0.1:8080/v1 `
  --api-key-env LOCAL_LLM_API_KEY `
  --model local-gguf `
  --data-dir .\data\chembl_3prop_2pos `
  --split test `
  --max-rows 20 `
  --few-shot-k 0 `
  --temperature 0 `
  --out .\reports\llm_edit_baselines\chembl_deepseek14b_zeroshot.jsonl
```

Evaluate:

```powershell
python .\scripts\run_llm_edit_baseline.py evaluate `
  --data-dir .\data\chembl_3prop_2pos `
  --split test `
  --max-rows 20 `
  --predictions .\reports\llm_edit_baselines\chembl_deepseek14b_zeroshot.jsonl `
  --out .\reports\llm_edit_baselines\chembl_deepseek14b_zeroshot_metrics.json
```

Run the three ADMET edit datasets:

```powershell
foreach ($d in @("chembl_3prop_2pos", "pubchem_3prop_2pos", "papyrus_3prop_2pos")) {
  python .\scripts\run_llm_edit_baseline.py run `
    --provider openai-compatible `
    --base-url http://127.0.0.1:8080/v1 `
    --api-key-env LOCAL_LLM_API_KEY `
    --model local-gguf `
    --data-dir ".\data\$d" `
    --split test `
    --max-rows 50 `
    --few-shot-k 2 `
    --temperature 0 `
    --out ".\reports\llm_edit_baselines\${d}_local_fewshot.jsonl"

  python .\scripts\run_llm_edit_baseline.py evaluate `
    --data-dir ".\data\$d" `
    --split test `
    --max-rows 50 `
    --predictions ".\reports\llm_edit_baselines\${d}_local_fewshot.jsonl" `
    --out ".\reports\llm_edit_baselines\${d}_local_fewshot_metrics.json"
}
```

## 5. Run Zero-Shot / Few-Shot On BindingDB

Seen-target test:

```powershell
python .\scripts\run_llm_bindingdb_baseline.py run `
  --provider openai-compatible `
  --base-url http://127.0.0.1:8080/v1 `
  --api-key-env LOCAL_LLM_API_KEY `
  --model local-gguf `
  --data-dir .\data\bindingdb_target_conditioned `
  --split test_seen_target `
  --max-rows 50 `
  --few-shot-k 2 `
  --temperature 0 `
  --out .\reports\llm_bindingdb_baselines\bindingdb_seen_local_fewshot.jsonl

python .\scripts\run_llm_bindingdb_baseline.py evaluate `
  --data-dir .\data\bindingdb_target_conditioned `
  --split test_seen_target `
  --max-rows 50 `
  --predictions .\reports\llm_bindingdb_baselines\bindingdb_seen_local_fewshot.jsonl `
  --out .\reports\llm_bindingdb_baselines\bindingdb_seen_local_fewshot_metrics.json
```

Unseen-target test:

```powershell
python .\scripts\run_llm_bindingdb_baseline.py run `
  --provider openai-compatible `
  --base-url http://127.0.0.1:8080/v1 `
  --api-key-env LOCAL_LLM_API_KEY `
  --model local-gguf `
  --data-dir .\data\bindingdb_target_conditioned `
  --split test_unseen_target `
  --max-rows 50 `
  --few-shot-k 2 `
  --temperature 0 `
  --out .\reports\llm_bindingdb_baselines\bindingdb_unseen_local_fewshot.jsonl

python .\scripts\run_llm_bindingdb_baseline.py evaluate `
  --data-dir .\data\bindingdb_target_conditioned `
  --split test_unseen_target `
  --max-rows 50 `
  --predictions .\reports\llm_bindingdb_baselines\bindingdb_unseen_local_fewshot.jsonl `
  --out .\reports\llm_bindingdb_baselines\bindingdb_unseen_local_fewshot_metrics.json
```

## 6. LoRA / QLoRA Fine-Tuning

Do not fine-tune GGUF files directly. GGUF is for inference. For LoRA, use the original Hugging Face `safetensors` model, for example:

```text
Qwen/Qwen2.5-7B-Instruct
Qwen/Qwen2.5-14B-Instruct
deepseek-ai/DeepSeek-R1-Distill-Qwen-14B
```

Practical recommendation:

| Hardware | Recommendation |
|---|---|
| 16GB system RAM, no CUDA GPU | Do not fine-tune locally. Run GGUF zero-shot/few-shot only. |
| 12-16GB VRAM | Try Qwen2.5-7B-Instruct with QLoRA, batch size 1, gradient accumulation. |
| 24GB+ VRAM | Qwen2.5-14B or DeepSeek-R1-Distill-Qwen-14B QLoRA becomes more realistic. |

Install training dependencies in a GPU-capable environment:

```powershell
python -m venv .venv_ft
.\.venv_ft\Scripts\activate
python -m pip install -U pip
pip install torch transformers datasets peft accelerate bitsandbytes
```

Train ADMET merged LoRA:

```powershell
python .\scripts\train_lora_sft.py `
  --base-model Qwen/Qwen2.5-7B-Instruct `
  --train-jsonl .\data\llm_sft\admet_edit_3prop_merged\train.jsonl `
  --val-jsonl .\data\llm_sft\admet_edit_3prop_merged\val.jsonl `
  --out-dir .\models\lora\qwen2.5-7b-posway-admet-edit `
  --load-in-4bit `
  --gradient-checkpointing `
  --bf16 `
  --max-seq-length 2048 `
  --num-train-epochs 1 `
  --per-device-train-batch-size 1 `
  --gradient-accumulation-steps 16 `
  --learning-rate 2e-4
```

Train BindingDB pilot LoRA:

```powershell
python .\scripts\train_lora_sft.py `
  --base-model Qwen/Qwen2.5-7B-Instruct `
  --train-jsonl .\data\llm_sft_bindingdb_20k\bindingdb_target_conditioned\train.jsonl `
  --val-jsonl .\data\llm_sft_bindingdb_20k\bindingdb_target_conditioned\val.jsonl `
  --out-dir .\models\lora\qwen2.5-7b-posway-bindingdb-20k `
  --load-in-4bit `
  --gradient-checkpointing `
  --bf16 `
  --max-seq-length 2048 `
  --num-train-epochs 1 `
  --per-device-train-batch-size 1 `
  --gradient-accumulation-steps 16 `
  --learning-rate 2e-4
```

If BF16 is not supported, replace `--bf16` with `--fp16`.

## 7. Reporting Baselines In The Paper

Keep these rows separate:

| Baseline | Training? | Dataset |
|---|---|---|
| GGUF zero-shot | No | ChEMBL, PubChem, Papyrus, BindingDB. |
| GGUF few-shot | No | Same model, in-context examples only. |
| Qwen LoRA ADMET edit | Yes | Merged ChEMBL + PubChem + Papyrus only. |
| Qwen LoRA BindingDB | Yes | BindingDB target-conditioned only. |
| Traditional ranker | Yes | Candidate ranking baseline from `BASELINE_GUIDE.md`. |

Do not report a mixed ADMET+BindingDB fine-tune as the first result. Only do multi-task training after single-task baselines are established.
