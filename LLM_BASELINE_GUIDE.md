# LLM Direct-Edit Baseline Guide

This guide is for the baseline that matches the LLM task directly:

```text
input SMILES + natural-language instruction -> exactly 2 edited SMILES
```

This is different from the traditional candidate ranker in `BASELINE_GUIDE.md`.

## Script

```text
scripts/run_llm_edit_baseline.py
```

The script supports three steps:

| Step | Command | Purpose |
|---|---|---|
| Export prompts | `export-prompts` | Inspect or archive the exact prompts without calling a model. |
| Run model | `run` | Call GPT/OpenAI-compatible APIs or a local HuggingFace/Qwen model. |
| Evaluate | `evaluate` | Parse model outputs and compute edit metrics. |

## Evaluation Metrics

The evaluator reports:

- valid SMILES rate
- duplicate prediction rate
- exact hit@1
- exact hit@2
- exact recall@2 against the two gold positives
- all-gold exact@2
- mean best pair Tanimoto to gold positives
- mean gold coverage Tanimoto
- local constraint pass rate using MW, LogP, and QED default thresholds

Exact-match metrics will be harsh. Tanimoto metrics are important because a valid edited analog may be chemically close to gold even when it is not exactly one of the stored positives.

## Export Prompts First

Use this before spending API budget:

```powershell
python .\scripts\run_llm_edit_baseline.py export-prompts `
  --data-dir .\data\chembl_3prop_2pos `
  --split test `
  --max-rows 20 `
  --few-shot-k 2 `
  --out .\reports\llm_edit_baselines\chembl_test_prompts.jsonl
```

## GPT / OpenAI-Compatible API Baseline

Set your API key:

```powershell
$env:OPENAI_API_KEY="YOUR_KEY"
```

Run zero-shot:

```powershell
python .\scripts\run_llm_edit_baseline.py run `
  --provider openai-compatible `
  --base-url https://api.openai.com/v1 `
  --api-key-env OPENAI_API_KEY `
  --model YOUR_GPT_MODEL `
  --data-dir .\data\chembl_3prop_2pos `
  --split test `
  --max-rows 50 `
  --few-shot-k 0 `
  --temperature 0 `
  --out .\reports\llm_edit_baselines\chembl_gpt_zeroshot.jsonl
```

Run few-shot:

```powershell
python .\scripts\run_llm_edit_baseline.py run `
  --provider openai-compatible `
  --base-url https://api.openai.com/v1 `
  --api-key-env OPENAI_API_KEY `
  --model YOUR_GPT_MODEL `
  --data-dir .\data\chembl_3prop_2pos `
  --split test `
  --max-rows 50 `
  --few-shot-k 3 `
  --temperature 0 `
  --out .\reports\llm_edit_baselines\chembl_gpt_fewshot.jsonl
```

Evaluate:

```powershell
python .\scripts\run_llm_edit_baseline.py evaluate `
  --data-dir .\data\chembl_3prop_2pos `
  --split test `
  --max-rows 50 `
  --predictions .\reports\llm_edit_baselines\chembl_gpt_fewshot.jsonl `
  --out .\reports\llm_edit_baselines\chembl_gpt_fewshot_metrics.json
```

## Qwen / Local HuggingFace Baseline

Use a local model path or HuggingFace model ID that is already available in your environment/cache:

```powershell
python .\scripts\run_llm_edit_baseline.py run `
  --provider hf-local `
  --model Qwen/Qwen2.5-7B-Instruct `
  --data-dir .\data\chembl_3prop_2pos `
  --split test `
  --max-rows 20 `
  --few-shot-k 2 `
  --temperature 0 `
  --max-tokens 256 `
  --out .\reports\llm_edit_baselines\chembl_qwen_fewshot.jsonl
```

If the model requires custom code:

```powershell
  --trust-remote-code
```

If GPU memory is tight, start with fewer rows and a smaller model. This script does not download or install models by itself.

## Qwen Fine-Tuning Data

The prompt export JSONL can be converted to SFT format. Each line already contains:

- `messages`
- `gold_smiles`
- `input_smiles`
- `instruction`

For SFT, use `messages` as the conversation and replace the final assistant answer with:

```json
{"edited_smiles":["gold_smiles_1","gold_smiles_2"]}
```

A fine-tuned Qwen baseline should be reported separately from zero-shot/few-shot Qwen:

| Baseline | Meaning |
|---|---|
| Qwen zero-shot | No task examples. |
| Qwen few-shot | In-context train examples only. |
| Qwen SFT/LoRA | Supervised fine-tuned on train split. |

## Paper-Framing Recommendation

For the first serious table, run:

1. GPT zero-shot direct edit.
2. GPT few-shot direct edit.
3. Qwen zero-shot direct edit.
4. Qwen few-shot direct edit.
5. Traditional supervised ranker from `BASELINE_GUIDE.md`.

Then the paper can fairly answer:

```text
Do general/chemistry-capable LLMs actually edit molecules under multi-property constraints,
or do traditional supervised candidate rankers remain stronger on small experimental datasets?
```
