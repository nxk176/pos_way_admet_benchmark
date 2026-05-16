# Traditional Baseline Guide

This guide is for paper-oriented baselines that match the current task framing:

```text
conditional molecule editing / candidate ranking
```

These baselines are not diffusion, not de novo generation, and not LLM-only prompting. They are supervised or classical candidate rankers using molecule graph features plus instruction/metadata features.

## Script

```text
scripts/train_traditional_baselines.py
```

The script has two modes:

| Mode | Dataset | Task |
|---|---|---|
| `edit-ranker` | ChEMBL / PubChem / Papyrus 3-property datasets | Given input molecule + instruction + candidate pool, rank the correct edited molecules above decoys. |
| `bindingdb-ranker` | BindingDB target-conditioned triplets | Given input ligand + target context + positive/negative candidates, rank the stronger ligand higher. |

## Features

The supervised ranker uses:

- Morgan fingerprint bits for input molecule.
- Morgan fingerprint bits for candidate molecule.
- XOR/difference fingerprint bits between input and candidate.
- Hashed instruction tokens.
- Hashed metadata such as source, endpoint, direction, target, and measurement type.
- Simple numeric molecule features: Tanimoto similarity, MW delta, LogP delta, QED delta, heavy atom delta, candidate MW/LogP, and evidence count when available.

The model is `sklearn.linear_model.SGDClassifier(loss="log_loss")`, so it is a lightweight classical supervised baseline.

Important fairness note: the BindingDB ranker does **not** use p-values or positive/negative deltas as input features. Those fields are labels/evaluation evidence only.

## Dependencies

Use an environment with:

- RDKit
- numpy
- scipy
- scikit-learn
- joblib

From this workspace, the likely environment is:

```powershell
..\myenv311\Scripts\python.exe
```

## ChEMBL Smoke Test

Run this first to check the pipeline quickly:

```powershell
..\myenv311\Scripts\python.exe .\scripts\train_traditional_baselines.py edit-ranker `
  --data-dir .\data\chembl_3prop_2pos `
  --source-name ChEMBL `
  --max-train-queries 500 `
  --max-eval-queries 100 `
  --decoys-per-query 20 `
  --max-iter 10
```

## ChEMBL Full Baseline

```powershell
..\myenv311\Scripts\python.exe .\scripts\train_traditional_baselines.py edit-ranker `
  --data-dir .\data\chembl_3prop_2pos `
  --source-name ChEMBL `
  --decoys-per-query 20 `
  --max-iter 30 `
  --save-model
```

## PubChem Sample Baseline

PubChem is larger; start sampled:

```powershell
..\myenv311\Scripts\python.exe .\scripts\train_traditional_baselines.py edit-ranker `
  --data-dir .\data\pubchem_3prop_2pos `
  --source-name PubChem `
  --max-train-queries 5000 `
  --max-eval-queries 1000 `
  --decoys-per-query 20 `
  --max-iter 20
```

## Papyrus Sample Baseline

```powershell
..\myenv311\Scripts\python.exe .\scripts\train_traditional_baselines.py edit-ranker `
  --data-dir .\data\papyrus_3prop_2pos `
  --source-name Papyrus `
  --max-train-queries 5000 `
  --max-eval-queries 1000 `
  --decoys-per-query 20 `
  --max-iter 20
```

## BindingDB Smoke Test

```powershell
..\myenv311\Scripts\python.exe .\scripts\train_traditional_baselines.py bindingdb-ranker `
  --data-dir .\data\bindingdb_target_conditioned `
  --max-train-rows 20000 `
  --max-eval-rows 5000 `
  --max-iter 15
```

## BindingDB Larger Run

```powershell
..\myenv311\Scripts\python.exe .\scripts\train_traditional_baselines.py bindingdb-ranker `
  --data-dir .\data\bindingdb_target_conditioned `
  --max-train-rows 0 `
  --max-eval-rows 0 `
  --max-iter 25 `
  --save-model
```

## Outputs

Each run writes a timestamped folder under:

```text
reports/traditional_baselines/
```

Files:

- `metrics.json`: main results.
- `config.json`: exact command configuration.
- `model.joblib`: saved only when `--save-model` is passed.

## Metrics

For 3-property edit datasets:

- `hit_at_1`
- `recall_at_2`
- `all_gold_at_2`
- `mrr`

For BindingDB:

- `pairwise_accuracy`
- `candidate_auc`
- `candidate_average_precision`

The report includes a similarity baseline and the supervised SGD ranker, so the first paper question is:

```text
Does a cheap traditional supervised ranker beat similarity/MMP-style ranking?
```

If yes, that supports the boss's hypothesis that traditional baselines may be strong enough and that expensive generative models are not automatically necessary.
