# Local LLM Models For Baseline Runs

This file tracks local/open-weight models intended for the LLM direct-edit baseline:

```text
input SMILES + instruction -> exactly 2 edited SMILES
```

## Current Local Files

Current GGUF files found directly under `models/`:

| Model | Local path | Role |
|---|---|---|
| DeepSeek-R1-Distill-Qwen-14B Q3_K_M | `models/DeepSeek-R1-Distill-Qwen-14B-Q3_K_M.gguf` | Main 14B local reasoning-distilled baseline; single GGUF file, ready to host. |
| Qwen2.5-14B-Instruct Q3_K_M | `models/qwen2.5-14b-instruct-q3_k_m-00001-of-00002.gguf` + `00002-of-00002` | Official Qwen 14B instruct baseline; both shards are present. Point llama.cpp to shard `00001-of-00002`. |
| Qwen2.5-7B-Instruct Q5_K_M | `models/qwen2.5-7b-instruct-q5_k_m-00001-of-00002.gguf` + `00002-of-00002` | Official Qwen 7B instruct baseline; both shards are present. Point llama.cpp to shard `00001-of-00002`. |

Important: `.gguf` files are for inference. For LoRA/QLoRA fine-tuning, use the original Hugging Face `safetensors` model, train an adapter, then optionally export/convert for local GGUF inference.

## Recommended Baselines For 16GB RAM

| Priority | Model | Quantization | Why |
|---:|---|---|---|
| 1 | Qwen2.5-7B-Instruct | Q4_K_M | Official Qwen GGUF, Apache-2.0, strong local instruct baseline. |
| 2 | DeepSeek-R1-Distill-Qwen-14B | Q3_K_M or Q4_K_M | Reasoning-distilled model; 14B is more credible than 7B if RAM allows. |
| 3 | Qwen2.5-14B-Instruct | Q3_K_M or Q4_K_M | Official Qwen 14B instruct; heavier, useful if 16GB RAM can handle short context. |

Practical note for 16GB RAM:

- 7B Q4_K_M should be comfortable.
- 14B Q4_K_M can be tight once context/KV cache is included.
- 14B Q3_K_M is safer for local CPU/RAM runs.
- Use short context first, e.g. 4096 tokens.

## Direct Links

Official Qwen 7B GGUF:

- Repo: https://huggingface.co/Qwen/Qwen2.5-7B-Instruct-GGUF
- Q4_K_M shard 1: https://huggingface.co/Qwen/Qwen2.5-7B-Instruct-GGUF/resolve/main/qwen2.5-7b-instruct-q4_k_m-00001-of-00002.gguf?download=true
- Q4_K_M shard 2: https://huggingface.co/Qwen/Qwen2.5-7B-Instruct-GGUF/resolve/main/qwen2.5-7b-instruct-q4_k_m-00002-of-00002.gguf?download=true

DeepSeek-R1-Distill-Qwen-14B GGUF:

- Base official model: https://huggingface.co/deepseek-ai/DeepSeek-R1-Distill-Qwen-14B
- GGUF repo: https://huggingface.co/bartowski/DeepSeek-R1-Distill-Qwen-14B-GGUF
- Safer 16GB option, Q3_K_M: https://huggingface.co/bartowski/DeepSeek-R1-Distill-Qwen-14B-GGUF/resolve/main/DeepSeek-R1-Distill-Qwen-14B-Q3_K_M.gguf?download=true
- Higher-quality but tighter, Q4_K_M: https://huggingface.co/bartowski/DeepSeek-R1-Distill-Qwen-14B-GGUF/resolve/main/DeepSeek-R1-Distill-Qwen-14B-Q4_K_M.gguf?download=true

Official Qwen 14B GGUF:

- Repo: https://huggingface.co/Qwen/Qwen2.5-14B-Instruct-GGUF
- Q4_K_M shard 1: https://huggingface.co/Qwen/Qwen2.5-14B-Instruct-GGUF/resolve/main/qwen2.5-14b-instruct-q4_k_m-00001-of-00003.gguf?download=true
- Q4_K_M shard 2: https://huggingface.co/Qwen/Qwen2.5-14B-Instruct-GGUF/resolve/main/qwen2.5-14b-instruct-q4_k_m-00002-of-00003.gguf?download=true
- Q4_K_M shard 3: https://huggingface.co/Qwen/Qwen2.5-14B-Instruct-GGUF/resolve/main/qwen2.5-14b-instruct-q4_k_m-00003-of-00003.gguf?download=true

## Resume Download Commands

DeepSeek-R1-Distill-Qwen-14B Q3_K_M, safer on 16GB:

```powershell
New-Item -ItemType Directory -Force .\models\gguf\DeepSeek-R1-Distill-Qwen-14B | Out-Null
curl.exe -L -C - `
  -o .\models\gguf\DeepSeek-R1-Distill-Qwen-14B\DeepSeek-R1-Distill-Qwen-14B-Q3_K_M.gguf `
  "https://huggingface.co/bartowski/DeepSeek-R1-Distill-Qwen-14B-GGUF/resolve/main/DeepSeek-R1-Distill-Qwen-14B-Q3_K_M.gguf?download=true"
```

DeepSeek-R1-Distill-Qwen-14B Q4_K_M, better but tighter:

```powershell
New-Item -ItemType Directory -Force .\models\gguf\DeepSeek-R1-Distill-Qwen-14B | Out-Null
curl.exe -L -C - `
  -o .\models\gguf\DeepSeek-R1-Distill-Qwen-14B\DeepSeek-R1-Distill-Qwen-14B-Q4_K_M.gguf `
  "https://huggingface.co/bartowski/DeepSeek-R1-Distill-Qwen-14B-GGUF/resolve/main/DeepSeek-R1-Distill-Qwen-14B-Q4_K_M.gguf?download=true"
```

Qwen2.5-14B-Instruct Q4_K_M, official Qwen but split into 3 shards:

```powershell
New-Item -ItemType Directory -Force .\models\gguf\Qwen2.5-14B-Instruct | Out-Null
curl.exe -L -C - -o .\models\gguf\Qwen2.5-14B-Instruct\qwen2.5-14b-instruct-q4_k_m-00001-of-00003.gguf "https://huggingface.co/Qwen/Qwen2.5-14B-Instruct-GGUF/resolve/main/qwen2.5-14b-instruct-q4_k_m-00001-of-00003.gguf?download=true"
curl.exe -L -C - -o .\models\gguf\Qwen2.5-14B-Instruct\qwen2.5-14b-instruct-q4_k_m-00002-of-00003.gguf "https://huggingface.co/Qwen/Qwen2.5-14B-Instruct-GGUF/resolve/main/qwen2.5-14b-instruct-q4_k_m-00002-of-00003.gguf?download=true"
curl.exe -L -C - -o .\models\gguf\Qwen2.5-14B-Instruct\qwen2.5-14b-instruct-q4_k_m-00003-of-00003.gguf "https://huggingface.co/Qwen/Qwen2.5-14B-Instruct-GGUF/resolve/main/qwen2.5-14b-instruct-q4_k_m-00003-of-00003.gguf?download=true"
```

## Serving With llama.cpp

If using llama.cpp, the easiest route is often to let `llama-server` pull from Hugging Face directly:

```powershell
llama-server -hf Qwen/Qwen2.5-7B-Instruct-GGUF:Q4_K_M -c 4096 --port 8080
```

For the 14B distill:

```powershell
llama-server -hf bartowski/DeepSeek-R1-Distill-Qwen-14B-GGUF:Q3_K_M -c 4096 --port 8080
```

Then run the local OpenAI-compatible baseline:

```powershell
$env:LOCAL_LLM_API_KEY="dummy"
python .\scripts\run_llm_edit_baseline.py run `
  --provider openai-compatible `
  --base-url http://127.0.0.1:8080/v1 `
  --api-key-env LOCAL_LLM_API_KEY `
  --model local-llm `
  --data-dir .\data\chembl_3prop_2pos `
  --split test `
  --max-rows 20 `
  --few-shot-k 2 `
  --temperature 0 `
  --out .\reports\llm_edit_baselines\chembl_local_llm_fewshot.jsonl
```
