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
    from rdkit.Chem import AllChem, Descriptors, QED
except Exception:  # noqa: BLE001
    Chem = None
    DataStructs = None
    RDLogger = None
    AllChem = None
    Descriptors = None
    QED = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate LoRA JSON predictions for ADMET edit or BindingDB SFT data.")
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--reference-jsonl",
        type=Path,
        action="append",
        default=[],
        help="Optional reference SFT JSONL file. Can be repeated. Used to recover metadata not stored in predictions.",
    )
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


def load_reference_rows(paths: list[Path]) -> dict[str, dict[str, Any]]:
    refs: dict[str, dict[str, Any]] = {}
    for path in paths:
        for row in read_jsonl(path):
            row_id = str(row.get("id") or "")
            if row_id:
                refs[row_id] = row
    return refs


def merge_reference_metadata(row: dict[str, Any], refs: dict[str, dict[str, Any]]) -> dict[str, Any]:
    ref = refs.get(str(row.get("id") or ""))
    if not ref:
        return row
    merged = dict(row)
    for key in [
        "split",
        "primary_endpoint",
        "primary_objective",
        "preserved_property",
        "local_constraints",
        "positive_answers",
        "source_positive_sample_ids",
    ]:
        value = merged.get(key)
        if value in (None, "", []):
            merged[key] = ref.get(key)
    return merged


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


def mol_descriptors(smiles: str) -> dict[str, float] | None:
    if Chem is None or Descriptors is None or QED is None:
        return None
    mol = Chem.MolFromSmiles(str(smiles or "").strip())
    if mol is None:
        return None
    return {
        "mw": float(Descriptors.MolWt(mol)),
        "logp": float(Descriptors.MolLogP(mol)),
        "qed": float(QED.qed(mol)),
    }


def parse_model_payload(text: str | None) -> dict[str, Any] | None:
    text = text or ""
    match = JSON_RE.search(text)
    if not match:
        start = text.find("{")
        if start < 0:
            return None
        candidate = text[start:]
    else:
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


def mean_or_none(values: list[float]) -> float | None:
    return round(sum(values) / len(values), 6) if values else None


def to_float(value: Any) -> float | None:
    try:
        if value in (None, ""):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


SUMMARY_KEYS = [
    "valid_smiles_rate",
    "exact_hit_any",
    "exact_recall",
    "mean_best_tanimoto",
    "rows_with_prediction_rate",
    "rows_with_two_or_more_predictions_rate",
    "rows_with_exactly_two_predictions_rate",
    "unchanged_input_rate",
    "rows_with_unchanged_input_rate",
    "duplicate_output_row_rate",
    "mean_num_predictions",
    "mean_num_valid_predictions",
    "mean_best_input_tanimoto",
    "mean_prediction_input_tanimoto",
    "descriptor_mw_delta_lte_100_rate",
    "descriptor_logp_delta_lte_2_rate",
    "descriptor_qed_drop_lte_0_1_rate",
    "descriptor_3check_pass_rate",
    "mean_abs_delta_mw",
    "mean_abs_delta_logp",
    "mean_delta_qed",
    "gold_primary_delta_mean",
    "exact_matched_primary_delta_mean",
]


def summarize_metric_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {"rows": len(rows)}
    for key in SUMMARY_KEYS:
        values = [float(row[key]) for row in rows if row.get(key) is not None]
        out[key] = mean(values)
    return out


def macro_average(groups: dict[str, dict[str, Any]]) -> dict[str, Any]:
    if not groups:
        return {}
    return {
        "groups": len(groups),
        **{
            key: mean([float(payload[key]) for payload in groups.values() if payload.get(key) is not None])
            for key in SUMMARY_KEYS
        },
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


def positive_answer_records(row: dict[str, Any]) -> list[dict[str, Any]]:
    values = row.get("positive_answers")
    if not isinstance(values, list):
        return []
    return [item for item in values if isinstance(item, dict)]


def gold_delta_values(row: dict[str, Any]) -> list[float]:
    values = []
    for item in positive_answer_records(row):
        value = to_float(item.get("delta_value"))
        if value is not None:
            values.append(value)
    return values


def exact_matched_deltas(row: dict[str, Any], pred_canon: list[str]) -> list[float]:
    pred_set = set(pred_canon)
    values = []
    for item in positive_answer_records(row):
        target = canonical(str(item.get("target_smiles_canon") or ""))
        value = to_float(item.get("delta_value"))
        if target and target in pred_set and value is not None:
            values.append(value)
    return values


def main() -> int:
    if RDLogger is not None:
        RDLogger.DisableLog("rdApp.warning")
        RDLogger.DisableLog("rdApp.error")
    args = parse_args()
    refs = load_reference_rows(args.reference_jsonl)
    rows = [merge_reference_metadata(row, refs) for row in read_jsonl(args.predictions)]
    counters: Counter[str] = Counter()
    summary_rows = []
    per_row = []
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)

    for row in rows:
        preds = predicted_smiles(row)
        gold = row.get("gold_smiles") or []
        if isinstance(gold, str):
            gold = [gold]
        input_smiles = str(row.get("input_smiles") or "").strip()
        input_canon = canonical(input_smiles)
        pred_canon = [canonical(item) for item in preds]
        pred_canon = [item for item in pred_canon if item]
        gold_canon = [canonical(item) for item in gold]
        gold_canon = [item for item in gold_canon if item]
        if not preds:
            counters["empty_prediction"] += 1
        if len(pred_canon) < len(preds):
            counters["invalid_prediction_smiles"] += len(preds) - len(pred_canon)
        valid_rate = len(pred_canon) / max(len(preds), 1)
        matches = set(pred_canon) & set(gold_canon)
        exact_any_value = 1.0 if matches else 0.0
        exact_recall_value = len(matches) / max(len(gold_canon), 1)
        sims = [
            sim
            for p in pred_canon
            for g in gold_canon
            for sim in [tanimoto(p, g)]
            if sim is not None
        ]
        input_sims = [
            sim
            for p in pred_canon
            for sim in [tanimoto(p, input_canon)]  # type: ignore[arg-type]
            if input_canon and sim is not None
        ]
        unchanged_count = sum(1 for item in pred_canon if input_canon and item == input_canon)
        if unchanged_count:
            counters["unchanged_input_prediction"] += unchanged_count
        duplicate_output = bool(pred_canon and len(pred_canon) != len(set(pred_canon)))
        if duplicate_output:
            counters["duplicate_output_row"] += 1

        input_desc = mol_descriptors(input_canon or "")
        desc_rows = []
        if input_desc is not None:
            for item in pred_canon:
                pred_desc = mol_descriptors(item)
                if pred_desc is None:
                    continue
                delta_mw = pred_desc["mw"] - input_desc["mw"]
                delta_logp = pred_desc["logp"] - input_desc["logp"]
                delta_qed = pred_desc["qed"] - input_desc["qed"]
                desc_rows.append(
                    {
                        "abs_delta_mw": abs(delta_mw),
                        "abs_delta_logp": abs(delta_logp),
                        "delta_qed": delta_qed,
                        "mw_pass": abs(delta_mw) <= 100.0,
                        "logp_pass": abs(delta_logp) <= 2.0,
                        "qed_pass": delta_qed >= -0.1,
                    }
                )

        gold_deltas = gold_delta_values(row)
        matched_deltas = exact_matched_deltas(row, pred_canon)
        if row.get("task") == "admet_3property_2positive_edit" and not gold_deltas:
            counters["missing_gold_primary_delta_metadata"] += 1

        row_metrics = {
            "rows_with_prediction_rate": 1.0 if preds else 0.0,
            "rows_with_two_or_more_predictions_rate": 1.0 if len(preds) >= 2 else 0.0,
            "rows_with_exactly_two_predictions_rate": 1.0 if len(preds) == 2 else 0.0,
            "mean_num_predictions": float(len(preds)),
            "mean_num_valid_predictions": float(len(pred_canon)),
            "valid_smiles_rate": valid_rate,
            "exact_any": exact_any_value,
            "exact_hit_any": exact_any_value,
            "exact_recall": exact_recall_value,
            "best_tanimoto": max(sims) if sims else None,
            "mean_best_tanimoto": max(sims) if sims else None,
            "unchanged_input_rate": unchanged_count / max(len(pred_canon), 1),
            "rows_with_unchanged_input_rate": 1.0 if unchanged_count else 0.0,
            "duplicate_output_row_rate": 1.0 if duplicate_output else 0.0,
            "mean_best_input_tanimoto": max(input_sims) if input_sims else None,
            "mean_prediction_input_tanimoto": mean_or_none(input_sims),
            "descriptor_mw_delta_lte_100_rate": mean([float(item["mw_pass"]) for item in desc_rows])
            if desc_rows
            else None,
            "descriptor_logp_delta_lte_2_rate": mean([float(item["logp_pass"]) for item in desc_rows])
            if desc_rows
            else None,
            "descriptor_qed_drop_lte_0_1_rate": mean([float(item["qed_pass"]) for item in desc_rows])
            if desc_rows
            else None,
            "descriptor_3check_pass_rate": mean(
                [
                    float(bool(item["mw_pass"] and item["logp_pass"] and item["qed_pass"]))
                    for item in desc_rows
                ]
            )
            if desc_rows
            else None,
            "mean_abs_delta_mw": mean_or_none([float(item["abs_delta_mw"]) for item in desc_rows]),
            "mean_abs_delta_logp": mean_or_none([float(item["abs_delta_logp"]) for item in desc_rows]),
            "mean_delta_qed": mean_or_none([float(item["delta_qed"]) for item in desc_rows]),
            "gold_primary_delta_mean": mean_or_none(gold_deltas),
            "exact_matched_primary_delta_mean": mean_or_none(matched_deltas),
        }
        summary_rows.append(row_metrics)
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
                "best_input_tanimoto": max(input_sims) if input_sims else None,
                "num_predictions": len(preds),
                "num_valid_predictions": len(pred_canon),
                "unchanged_input_predictions": unchanged_count,
                "duplicate_output": duplicate_output,
                "descriptor_3check_pass_rate": row_metrics["descriptor_3check_pass_rate"],
                "gold_primary_deltas": gold_deltas,
                "exact_matched_primary_deltas": matched_deltas,
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
        "reference_jsonl": [str(path) for path in args.reference_jsonl],
        "rdkit_available": Chem is not None,
        "rows": len(rows),
        **{key: summarize_metric_rows(summary_rows).get(key) for key in SUMMARY_KEYS},
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
