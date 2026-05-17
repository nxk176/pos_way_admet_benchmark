from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import AllChem
except Exception:  # noqa: BLE001
    Chem = None
    DataStructs = None
    RDLogger = None
    AllChem = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate LoRA JSON predictions for ADMET edit or BindingDB SFT data.")
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    return parser.parse_args()


JSON_RE = re.compile(r"\{.*\}", re.DOTALL)
FIELD_RE = re.compile(r'"(?P<key>stronger_smiles|hard_negative_smiles|smiles)"\s*:\s*"(?P<value>[^"]*)"')
ARRAY_FIELD_RE = re.compile(r'"(?P<key>edited_smiles|molecules|smiles)"\s*:\s*\[(?P<value>.*?)\]', re.DOTALL)
ARRAY_ITEM_RE = re.compile(r'"([^"]*)"')


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def canonical(smiles: str) -> str | None:
    smiles = str(smiles or "").strip()
    if not smiles:
        return None
    if Chem is None:
        return smiles
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def parse_model_payload(text: str | None) -> dict[str, Any] | None:
    match = JSON_RE.search(text or "")
    if not match:
        return None
    candidate = match.group(0)
    try:
        payload = json.loads(candidate)
        return payload if isinstance(payload, dict) else None
    except json.JSONDecodeError:
        pass

    repaired: dict[str, Any] = {}
    for field_match in FIELD_RE.finditer(candidate):
        repaired[field_match.group("key")] = field_match.group("value").strip()
    for field_match in ARRAY_FIELD_RE.finditer(candidate):
        values = [item.strip() for item in ARRAY_ITEM_RE.findall(field_match.group("value")) if item.strip()]
        if values:
            repaired[field_match.group("key")] = values
    return repaired or None


def tanimoto(a: str, b: str) -> float | None:
    if Chem is None or AllChem is None or DataStructs is None:
        return None
    ma = Chem.MolFromSmiles(a)
    mb = Chem.MolFromSmiles(b)
    if ma is None or mb is None:
        return None
    fpa = AllChem.GetMorganFingerprintAsBitVect(ma, 2, nBits=2048)
    fpb = AllChem.GetMorganFingerprintAsBitVect(mb, 2, nBits=2048)
    return float(DataStructs.FingerprintSimilarity(fpa, fpb))


def mean(values: list[float]) -> float:
    return round(sum(values) / len(values), 6) if values else 0.0


def summarize_metric_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "rows": len(rows),
        "valid_smiles_rate": mean([float(row["valid_smiles_rate"]) for row in rows]),
        "exact_hit_any": mean([float(row["exact_any"]) for row in rows]),
        "exact_recall": mean([float(row["exact_recall"]) for row in rows]),
        "mean_best_tanimoto": mean(
            [float(row["best_tanimoto"]) for row in rows if row.get("best_tanimoto") is not None]
        ),
    }


def macro_average(groups: dict[str, dict[str, Any]]) -> dict[str, Any]:
    if not groups:
        return {}
    keys = ["valid_smiles_rate", "exact_hit_any", "exact_recall", "mean_best_tanimoto"]
    return {
        "groups": len(groups),
        **{key: mean([float(payload[key]) for payload in groups.values()]) for key in keys},
    }


def predicted_smiles(row: dict[str, Any]) -> list[str]:
    payload = row.get("parsed_json")
    if not isinstance(payload, dict):
        payload = parse_model_payload(row.get("raw_output"))
    if not isinstance(payload, dict):
        return []
    if row.get("task") == "bindingdb_target_conditioned_triplet":
        value = payload.get("stronger_smiles") or payload.get("edited_smiles") or payload.get("smiles")
        return [value] if isinstance(value, str) and value.strip() else []
    values = payload.get("edited_smiles") or payload.get("smiles") or payload.get("molecules")
    if isinstance(values, str):
        values = [values]
    if not isinstance(values, list):
        return []
    out = []
    for item in values:
        if isinstance(item, dict):
            item = item.get("smiles")
        if isinstance(item, str) and item.strip():
            out.append(item.strip())
    return out


def main() -> int:
    if RDLogger is not None:
        RDLogger.DisableLog("rdApp.warning")
        RDLogger.DisableLog("rdApp.error")
    args = parse_args()
    rows = read_jsonl(args.predictions)
    counters: Counter[str] = Counter()
    valid_rates = []
    exact_any = []
    exact_recall = []
    best_tanimoto = []
    per_row = []
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)

    for row in rows:
        preds = predicted_smiles(row)
        gold = row.get("gold_smiles") or []
        if isinstance(gold, str):
            gold = [gold]
        pred_canon = [canonical(item) for item in preds]
        pred_canon = [item for item in pred_canon if item]
        gold_canon = [canonical(item) for item in gold]
        gold_canon = [item for item in gold_canon if item]
        if not preds:
            counters["empty_prediction"] += 1
        valid_rate = len(pred_canon) / max(len(preds), 1)
        valid_rates.append(valid_rate)
        matches = set(pred_canon) & set(gold_canon)
        exact_any_value = 1.0 if matches else 0.0
        exact_recall_value = len(matches) / max(len(gold_canon), 1)
        exact_any.append(exact_any_value)
        exact_recall.append(exact_recall_value)
        sims = [
            sim
            for p in pred_canon
            for g in gold_canon
            for sim in [tanimoto(p, g)]
            if sim is not None
        ]
        if sims:
            best_tanimoto.append(max(sims))
        row_metrics = {
            "valid_smiles_rate": valid_rate,
            "exact_any": exact_any_value,
            "exact_recall": exact_recall_value,
            "best_tanimoto": max(sims) if sims else None,
        }
        source_dataset = str(row.get("source_dataset") or "unknown")
        task = str(row.get("task") or "unknown")
        grouped[f"source:{source_dataset}"].append(row_metrics)
        grouped[f"task:{task}"].append(row_metrics)
        per_row.append(
            {
                "id": row.get("id"),
                "task": row.get("task"),
                "source_dataset": row.get("source_dataset"),
                "predicted_smiles": preds,
                "gold_smiles": gold,
                "exact_any": bool(matches),
                "best_tanimoto": max(sims) if sims else None,
            }
        )

    by_source = {
        key.replace("source:", "", 1): summarize_metric_rows(value)
        for key, value in sorted(grouped.items())
        if key.startswith("source:")
    }
    by_task = {
        key.replace("task:", "", 1): summarize_metric_rows(value)
        for key, value in sorted(grouped.items())
        if key.startswith("task:")
    }

    metrics = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "predictions": str(args.predictions),
        "rdkit_available": Chem is not None,
        "rows": len(rows),
        "valid_smiles_rate": mean(valid_rates),
        "exact_hit_any": mean(exact_any),
        "exact_recall": mean(exact_recall),
        "mean_best_tanimoto": mean(best_tanimoto),
        "by_source_dataset": by_source,
        "source_macro_average": macro_average(by_source),
        "by_task": by_task,
        "counters": dict(sorted(counters.items())),
        "per_row": per_row,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(metrics, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in metrics.items() if k != "per_row"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
