from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import sqlite3
import sys
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from statistics import median
from typing import Any, Iterable

try:
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
    from rdkit.Chem.Scaffolds import MurckoScaffold
except Exception as exc:  # noqa: BLE001
    raise SystemExit(
        "ERROR: this script requires RDKit. Run with .\\myenv311\\Scripts\\python.exe"
    ) from exc

RDLogger.DisableLog("rdApp.warning")
RDLogger.DisableLog("rdApp.error")


MEASUREMENTS = {
    "Ki (nM)": ("Ki", "binding_affinity", "pKi"),
    "Kd (nM)": ("Kd", "binding_affinity", "pKd"),
    "IC50 (nM)": ("IC50", "activity_potency", "pIC50"),
    "EC50 (nM)": ("EC50", "activity_potency", "pEC50"),
}

OBSERVATION_FIELDS = [
    "observation_uid",
    "source",
    "source_version",
    "source_record_id",
    "ligand_id",
    "smiles_raw",
    "inchikey",
    "connectivity_key",
    "target_id",
    "target_name",
    "target_organism",
    "target_sequence_hash",
    "target_sequence_length",
    "measurement_type",
    "measurement_group",
    "relation",
    "value_nm",
    "p_value",
    "unit_raw",
    "unit_canonical",
    "condition_bucket",
    "ph",
    "temp_c",
    "curation_source",
    "doi",
    "pmid",
    "pubchem_aid",
    "patent_number",
    "activity_direction",
    "modulation_label",
    "modulation_numeric",
]

AGGREGATED_FIELDS = [
    "aggregate_uid",
    "source",
    "source_version",
    "ligand_id",
    "smiles_raw",
    "inchikey",
    "connectivity_key",
    "target_id",
    "target_name",
    "target_organism",
    "target_sequence_hash",
    "target_sequence_length",
    "measurement_type",
    "measurement_group",
    "p_value_median",
    "p_value_min",
    "p_value_max",
    "evidence_count",
    "pmid_count",
    "doi_count",
    "assay_count",
    "unit_canonical",
    "activity_direction",
    "modulation_label",
    "modulation_numeric",
]

RANK_READY_FIELDS = AGGREGATED_FIELDS + [
    "canonical_smiles",
    "rdkit_inchikey",
    "mw",
    "logp",
    "tpsa",
    "hba",
    "hbd",
    "rotatable_bonds",
    "heavy_atoms",
    "qed",
    "murcko_scaffold",
    "brics_fragments",
]

TRIPLET_FIELDS = [
    "sample_id",
    "source",
    "source_version",
    "task_type",
    "instruction",
    "target_id",
    "target_name",
    "target_organism",
    "target_sequence_hash",
    "target_sequence",
    "measurement_type",
    "measurement_group",
    "input_smiles",
    "input_inchikey",
    "input_connectivity_key",
    "input_p_value",
    "positive_smiles",
    "positive_inchikey",
    "positive_connectivity_key",
    "positive_p_value",
    "positive_delta",
    "negative_smiles",
    "negative_inchikey",
    "negative_connectivity_key",
    "negative_p_value",
    "negative_delta",
    "input_positive_tanimoto",
    "input_negative_tanimoto",
    "positive_negative_tanimoto",
    "shared_core",
    "shared_core_heavy_atoms",
    "input_mw",
    "positive_mw",
    "negative_mw",
    "input_logp",
    "positive_logp",
    "negative_logp",
    "input_qed",
    "positive_qed",
    "negative_qed",
    "input_evidence_count",
    "positive_evidence_count",
    "negative_evidence_count",
]

PAIR_FIELDS = [
    "pair_id",
    "source",
    "source_version",
    "target_id",
    "target_name",
    "target_sequence_hash",
    "measurement_type",
    "measurement_group",
    "input_smiles",
    "input_connectivity_key",
    "input_p_value",
    "positive_smiles",
    "positive_connectivity_key",
    "positive_p_value",
    "positive_delta",
    "input_positive_tanimoto",
    "shared_core",
    "input_evidence_count",
    "positive_evidence_count",
]

VALUE_RE = re.compile(r"([<>]=?|=)?\s*([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a target-conditioned BindingDB ranking/triplet dataset from a full BindingDB TSV snapshot."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("pos_way_admet_benchmark/raw/public/bindingdb/BindingDB_All_202605_tsv.zip"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/bindingdb_target_conditioned"),
    )
    parser.add_argument("--source-version", default="202605_full")
    parser.add_argument("--max-raw-rows", type=int, default=0, help="0 means full snapshot.")
    parser.add_argument("--min-positive-delta", type=float, default=0.5)
    parser.add_argument("--min-negative-delta", type=float, default=0.2)
    parser.add_argument("--min-input-candidate-tanimoto", type=float, default=0.35)
    parser.add_argument("--min-positive-negative-tanimoto", type=float, default=0.30)
    parser.add_argument("--min-scaffold-group-size", type=int, default=3)
    parser.add_argument("--max-candidate-scan", type=int, default=250)
    parser.add_argument("--max-queries-per-target-metric", type=int, default=0, help="0 means no cap.")
    parser.add_argument("--compute-brics-fragments", action="store_true")
    parser.add_argument("--compute-rdkit-inchikey", action="store_true")
    parser.add_argument("--compute-physchem-descriptors", action="store_true")
    parser.add_argument("--seen-val-bucket", type=int, default=10)
    parser.add_argument("--seen-test-bucket", type=int, default=10)
    parser.add_argument("--unseen-target-bucket", type=int, default=10)
    parser.add_argument("--keep-sqlite", action="store_true")
    return parser.parse_args()


def safe_field_size_limit() -> None:
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit //= 10


def clean_text(value: Any) -> str:
    return str(value or "").replace("\r", " ").replace("\n", " ").strip()


def slug(value: str, max_len: int = 80) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")
    return (cleaned or "unknown")[:max_len]


def stable_hash(value: str) -> int:
    digest = hashlib.sha1(value.encode("utf-8", errors="ignore")).hexdigest()
    return int(digest[:12], 16)


def short_hash(value: str, length: int = 16) -> str:
    return hashlib.sha1(value.encode("utf-8", errors="ignore")).hexdigest()[:length]


def connectivity_key(inchikey: str, smiles: str = "") -> str:
    inchikey = clean_text(inchikey)
    if inchikey:
        return inchikey.split("-", 1)[0]
    return f"SMILES_{short_hash(smiles, 14)}" if smiles else ""


def read_bindingdb_rows(path: Path) -> Iterable[dict[str, str]]:
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as zf:
            entries = [entry for entry in zf.infolist() if not entry.is_dir()]
            if not entries:
                raise SystemExit(f"ERROR: zip has no files: {path}")
            with zf.open(entries[0], "r") as raw:
                text = (line.decode("utf-8", errors="replace") for line in raw)
                yield from csv.DictReader(text, delimiter="\t")
    else:
        with path.open("r", encoding="utf-8", newline="") as handle:
            yield from csv.DictReader(handle, delimiter="\t")


def parse_nm_value(value: str) -> tuple[str, float, float] | None:
    text = clean_text(value).replace(",", "")
    if not text:
        return None
    match = VALUE_RE.search(text)
    if not match:
        return None
    relation = match.group(1) or "="
    nm_value = float(match.group(2))
    if not math.isfinite(nm_value) or nm_value <= 0.0:
        return None
    p_value = 9.0 - math.log10(nm_value)
    if not math.isfinite(p_value):
        return None
    return relation, nm_value, p_value


def target_from_row(row: dict[str, str]) -> dict[str, str]:
    num_chains_raw = clean_text(row.get("Number of Protein Chains in Target (>1 implies a multichain complex)"))
    try:
        num_chains = max(1, min(12, int(float(num_chains_raw)))) if num_chains_raw else 1
    except ValueError:
        num_chains = 1

    ids: list[str] = []
    sequences: list[str] = []
    names: list[str] = []
    for idx in range(1, num_chains + 1):
        seq = clean_text(row.get(f"BindingDB Target Chain Sequence {idx}"))
        sp_id = clean_text(row.get(f"UniProt (SwissProt) Primary ID of Target Chain {idx}"))
        trembl_id = clean_text(row.get(f"UniProt (TrEMBL) Primary ID of Target Chain {idx}"))
        rec_name = clean_text(row.get(f"UniProt (SwissProt) Recommended Name of Target Chain {idx}"))
        sub_name = clean_text(row.get(f"UniProt (TrEMBL) Submitted Name of Target Chain {idx}"))
        chain_id = sp_id or trembl_id or (f"SEQ_{short_hash(seq, 12)}" if seq else "")
        if chain_id:
            ids.append(chain_id)
        if seq:
            sequences.append(seq)
        if rec_name or sub_name:
            names.append(rec_name or sub_name)

    target_name = clean_text(row.get("Target Name")) or "; ".join(names)
    target_id = "+".join(ids) if ids else f"TARGET_{slug(target_name, 60)}"
    target_sequence = "|".join(sequences)
    sequence_hash = short_hash(target_sequence, 16) if target_sequence else ""
    return {
        "target_id": target_id,
        "target_name": target_name,
        "target_organism": clean_text(row.get("Target Source Organism According to Curator or DataSource")),
        "target_sequence": target_sequence,
        "target_sequence_hash": sequence_hash,
        "target_sequence_length": str(sum(len(seq) for seq in sequences)) if sequences else "",
    }


def ligand_id_from_row(row: dict[str, str]) -> str:
    monomer = clean_text(row.get("BindingDB MonomerID"))
    reactant_set = clean_text(row.get("BindingDB Reactant_set_id"))
    if monomer:
        return f"BDBM{monomer}"
    if reactant_set:
        return f"BDBR{reactant_set}"
    return ""


def init_db(db_path: Path) -> sqlite3.Connection:
    if db_path.exists():
        db_path.unlink()
    conn = sqlite3.connect(str(db_path))
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA temp_store=FILE")
    conn.execute(
        """
        CREATE TABLE exact_obs (
          connectivity_key TEXT,
          ligand_id TEXT,
          smiles_raw TEXT,
          inchikey TEXT,
          target_id TEXT,
          target_name TEXT,
          target_organism TEXT,
          target_sequence_hash TEXT,
          target_sequence_length TEXT,
          measurement_type TEXT,
          measurement_group TEXT,
          p_value REAL,
          unit_canonical TEXT,
          pmid TEXT,
          doi TEXT,
          pubchem_aid TEXT,
          activity_direction TEXT,
          modulation_label TEXT,
          modulation_numeric REAL
        )
        """
    )
    conn.execute(
        """
        CREATE TABLE targets (
          target_id TEXT PRIMARY KEY,
          target_name TEXT,
          target_organism TEXT,
          target_sequence_hash TEXT,
          target_sequence_length TEXT,
          target_sequence TEXT,
          exact_observation_count INTEGER DEFAULT 0
        )
        """
    )
    conn.execute(
        """
        CREATE TABLE aggregated (
          aggregate_uid TEXT PRIMARY KEY,
          source TEXT,
          source_version TEXT,
          ligand_id TEXT,
          smiles_raw TEXT,
          inchikey TEXT,
          connectivity_key TEXT,
          target_id TEXT,
          target_name TEXT,
          target_organism TEXT,
          target_sequence_hash TEXT,
          target_sequence_length TEXT,
          measurement_type TEXT,
          measurement_group TEXT,
          p_value_median REAL,
          p_value_min REAL,
          p_value_max REAL,
          evidence_count INTEGER,
          pmid_count INTEGER,
          doi_count INTEGER,
          assay_count INTEGER,
          unit_canonical TEXT,
          activity_direction TEXT,
          modulation_label TEXT,
          modulation_numeric REAL
        )
        """
    )
    conn.execute(
        """
        CREATE TABLE rank_ready (
          aggregate_uid TEXT PRIMARY KEY,
          source TEXT,
          source_version TEXT,
          ligand_id TEXT,
          smiles_raw TEXT,
          inchikey TEXT,
          connectivity_key TEXT,
          target_id TEXT,
          target_name TEXT,
          target_organism TEXT,
          target_sequence_hash TEXT,
          target_sequence_length TEXT,
          measurement_type TEXT,
          measurement_group TEXT,
          p_value_median REAL,
          p_value_min REAL,
          p_value_max REAL,
          evidence_count INTEGER,
          pmid_count INTEGER,
          doi_count INTEGER,
          assay_count INTEGER,
          unit_canonical TEXT,
          activity_direction TEXT,
          modulation_label TEXT,
          modulation_numeric REAL,
          canonical_smiles TEXT,
          rdkit_inchikey TEXT,
          mw REAL,
          logp REAL,
          tpsa REAL,
          hba INTEGER,
          hbd INTEGER,
          rotatable_bonds INTEGER,
          heavy_atoms INTEGER,
          qed REAL,
          murcko_scaffold TEXT,
          brics_fragments TEXT
        )
        """
    )
    return conn


def normalize_to_observations(args: argparse.Namespace, conn: sqlite3.Connection) -> dict[str, Any]:
    observations_path = args.out_dir / "bindingdb_observations.csv"
    counters: Counter[str] = Counter()
    exact_batch: list[tuple[Any, ...]] = []
    target_batch: dict[str, tuple[Any, ...]] = {}
    target_counts: Counter[str] = Counter()
    raw_rows = 0
    observation_rows = 0
    exact_rows = 0

    with observations_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OBSERVATION_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for row in read_bindingdb_rows(args.input):
            if args.max_raw_rows and raw_rows >= args.max_raw_rows:
                break
            raw_rows += 1
            smiles = clean_text(row.get("Ligand SMILES"))
            if not smiles:
                counters["skipped_missing_smiles"] += 1
                continue
            inchikey = clean_text(row.get("Ligand InChI Key"))
            conn_key = connectivity_key(inchikey, smiles)
            if not conn_key:
                counters["skipped_missing_connectivity"] += 1
                continue

            target = target_from_row(row)
            target_id = target["target_id"]
            target_batch[target_id] = (
                target_id,
                target["target_name"],
                target["target_organism"],
                target["target_sequence_hash"],
                target["target_sequence_length"],
                target["target_sequence"],
            )
            ligand_id = ligand_id_from_row(row)
            reactant_set_id = clean_text(row.get("BindingDB Reactant_set_id"))
            ph = clean_text(row.get("pH"))
            temp_c = clean_text(row.get("Temp (C)"))
            curation_source = clean_text(row.get("Curation/DataSource"))
            doi = clean_text(row.get("Article DOI")) or clean_text(row.get("BindingDB Entry DOI"))
            pmid = clean_text(row.get("PMID"))
            pubchem_aid = clean_text(row.get("PubChem AID"))
            patent_number = clean_text(row.get("Patent Number"))

            # BindingDB full TSV does not expose a reliable agonist/antagonist column.
            # Keep the boss-requested numeric coding but mark unknown as 0.5.
            modulation_label = "unknown"
            modulation_numeric = 0.5

            for column, (measurement_type, measurement_group, unit_canonical) in MEASUREMENTS.items():
                parsed = parse_nm_value(row.get(column, ""))
                if parsed is None:
                    continue
                relation, nm_value, p_value = parsed
                condition_bucket = "|".join(
                    [
                        target_id,
                        measurement_type,
                        f"pmid:{pmid or 'unknown'}",
                        f"doi:{doi or 'unknown'}",
                        f"aid:{pubchem_aid or 'unknown'}",
                    ]
                )
                obs_uid = f"bindingdb_{args.source_version}_{reactant_set_id or raw_rows}_{measurement_type}"
                out = {
                    "observation_uid": obs_uid,
                    "source": "BindingDB",
                    "source_version": args.source_version,
                    "source_record_id": reactant_set_id,
                    "ligand_id": ligand_id,
                    "smiles_raw": smiles,
                    "inchikey": inchikey,
                    "connectivity_key": conn_key,
                    "target_id": target_id,
                    "target_name": target["target_name"],
                    "target_organism": target["target_organism"],
                    "target_sequence_hash": target["target_sequence_hash"],
                    "target_sequence_length": target["target_sequence_length"],
                    "measurement_type": measurement_type,
                    "measurement_group": measurement_group,
                    "relation": relation,
                    "value_nm": round(nm_value, 8),
                    "p_value": round(p_value, 6),
                    "unit_raw": "nM",
                    "unit_canonical": unit_canonical,
                    "condition_bucket": condition_bucket,
                    "ph": ph,
                    "temp_c": temp_c,
                    "curation_source": curation_source,
                    "doi": doi,
                    "pmid": pmid,
                    "pubchem_aid": pubchem_aid,
                    "patent_number": patent_number,
                    "activity_direction": "higher_p_value_is_stronger",
                    "modulation_label": modulation_label,
                    "modulation_numeric": modulation_numeric,
                }
                writer.writerow(out)
                observation_rows += 1
                counters[f"{measurement_type}:parsed_observations"] += 1

                if relation != "=":
                    counters[f"{measurement_type}:censored_not_used_for_ranking"] += 1
                    continue

                exact_rows += 1
                target_counts[target_id] += 1
                exact_batch.append(
                    (
                        conn_key,
                        ligand_id,
                        smiles,
                        inchikey,
                        target_id,
                        target["target_name"],
                        target["target_organism"],
                        target["target_sequence_hash"],
                        target["target_sequence_length"],
                        measurement_type,
                        measurement_group,
                        float(p_value),
                        unit_canonical,
                        pmid,
                        doi,
                        pubchem_aid,
                        "higher_p_value_is_stronger",
                        modulation_label,
                        modulation_numeric,
                    )
                )
                if len(exact_batch) >= 5000:
                    conn.executemany(
                        "INSERT INTO exact_obs VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", exact_batch
                    )
                    exact_batch.clear()

            if raw_rows % 100000 == 0:
                print(f"[normalize] raw_rows={raw_rows:,} observations={observation_rows:,} exact={exact_rows:,}", flush=True)

    if exact_batch:
        conn.executemany("INSERT INTO exact_obs VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", exact_batch)
    conn.executemany(
        """
        INSERT OR IGNORE INTO targets (
          target_id,target_name,target_organism,target_sequence_hash,target_sequence_length,target_sequence,exact_observation_count
        ) VALUES (?,?,?,?,?,?,0)
        """,
        list(target_batch.values()),
    )
    for target_id, count in target_counts.items():
        conn.execute(
            "UPDATE targets SET exact_observation_count = exact_observation_count + ? WHERE target_id = ?",
            (int(count), target_id),
        )
    conn.commit()
    conn.execute("CREATE INDEX IF NOT EXISTS idx_exact_group ON exact_obs(connectivity_key,target_id,measurement_type)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_exact_target_metric ON exact_obs(target_id,measurement_type)")
    conn.commit()

    return {
        "raw_rows": raw_rows,
        "parsed_observations": observation_rows,
        "exact_observations_for_ranking": exact_rows,
        "counters": dict(sorted(counters.items())),
        "output": str(observations_path),
    }


def write_target_summary(args: argparse.Namespace, conn: sqlite3.Connection) -> dict[str, Any]:
    path = args.out_dir / "bindingdb_target_summary.csv"
    fields = [
        "target_id",
        "target_name",
        "target_organism",
        "target_sequence_hash",
        "target_sequence_length",
        "exact_observation_count",
        "target_sequence",
    ]
    rows = conn.execute(
        """
        SELECT target_id,target_name,target_organism,target_sequence_hash,target_sequence_length,
               exact_observation_count,target_sequence
        FROM targets
        ORDER BY exact_observation_count DESC, target_id
        """
    )
    count = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        for row in rows:
            writer.writerow(row)
            count += 1
    return {"targets": count, "output": str(path)}


def aggregate_exact_observations(args: argparse.Namespace, conn: sqlite3.Connection) -> dict[str, Any]:
    path = args.out_dir / "bindingdb_aggregated_observations.csv"
    insert_batch: list[tuple[Any, ...]] = []
    aggregate_rows = 0
    counters: Counter[str] = Counter()

    def flush_group(group_key: tuple[str, str, str] | None, records: list[sqlite3.Row], writer: csv.DictWriter) -> None:
        nonlocal aggregate_rows, insert_batch
        if not group_key or not records:
            return
        p_values = [float(r["p_value"]) for r in records]
        first = records[0]
        pmids = {clean_text(r["pmid"]) for r in records if clean_text(r["pmid"])}
        dois = {clean_text(r["doi"]) for r in records if clean_text(r["doi"])}
        assays = {clean_text(r["pubchem_aid"]) for r in records if clean_text(r["pubchem_aid"])}
        aggregate_uid = (
            f"bindingdbagg_{args.source_version}_"
            f"{short_hash(first['connectivity_key'] + '|' + first['target_id'] + '|' + first['measurement_type'], 18)}"
        )
        out = {
            "aggregate_uid": aggregate_uid,
            "source": "BindingDB",
            "source_version": args.source_version,
            "ligand_id": first["ligand_id"],
            "smiles_raw": first["smiles_raw"],
            "inchikey": first["inchikey"],
            "connectivity_key": first["connectivity_key"],
            "target_id": first["target_id"],
            "target_name": first["target_name"],
            "target_organism": first["target_organism"],
            "target_sequence_hash": first["target_sequence_hash"],
            "target_sequence_length": first["target_sequence_length"],
            "measurement_type": first["measurement_type"],
            "measurement_group": first["measurement_group"],
            "p_value_median": round(float(median(p_values)), 6),
            "p_value_min": round(min(p_values), 6),
            "p_value_max": round(max(p_values), 6),
            "evidence_count": len(records),
            "pmid_count": len(pmids),
            "doi_count": len(dois),
            "assay_count": len(assays),
            "unit_canonical": first["unit_canonical"],
            "activity_direction": first["activity_direction"],
            "modulation_label": first["modulation_label"],
            "modulation_numeric": first["modulation_numeric"],
        }
        writer.writerow(out)
        insert_batch.append(tuple(out[field] for field in AGGREGATED_FIELDS))
        if len(insert_batch) >= 5000:
            conn.executemany(f"INSERT INTO aggregated VALUES ({','.join(['?'] * len(AGGREGATED_FIELDS))})", insert_batch)
            insert_batch.clear()
        aggregate_rows += 1
        counters[f"{first['measurement_type']}:aggregated"] += 1

    conn.row_factory = sqlite3.Row
    query = """
        SELECT *
        FROM exact_obs
        ORDER BY connectivity_key, target_id, measurement_type
    """
    current_key: tuple[str, str, str] | None = None
    current_records: list[sqlite3.Row] = []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=AGGREGATED_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for row in conn.execute(query):
            key = (row["connectivity_key"], row["target_id"], row["measurement_type"])
            if current_key is not None and key != current_key:
                flush_group(current_key, current_records, writer)
                current_records = []
                if aggregate_rows and aggregate_rows % 100000 == 0:
                    print(f"[aggregate] aggregate_rows={aggregate_rows:,}", flush=True)
            current_key = key
            current_records.append(row)
        flush_group(current_key, current_records, writer)

    if insert_batch:
        conn.executemany(f"INSERT INTO aggregated VALUES ({','.join(['?'] * len(AGGREGATED_FIELDS))})", insert_batch)
    conn.commit()
    conn.execute("CREATE INDEX IF NOT EXISTS idx_agg_target_metric ON aggregated(target_id,measurement_type)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_agg_connectivity ON aggregated(connectivity_key)")
    conn.commit()
    return {"aggregated_rows": aggregate_rows, "counters": dict(sorted(counters.items())), "output": str(path)}


def descriptor_record(
    smiles: str,
    compute_brics_fragments: bool,
    compute_rdkit_inchikey: bool,
    compute_physchem_descriptors: bool,
) -> dict[str, Any] | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold = (
            Chem.MolToSmiles(scaffold_mol, canonical=True, isomericSmiles=False)
            if scaffold_mol is not None and scaffold_mol.GetNumHeavyAtoms() > 0
            else ""
        )
    except Exception:  # noqa: BLE001
        scaffold = ""
    brics_fragments = ""
    if compute_brics_fragments:
        try:
            from rdkit.Chem import BRICS

            brics_fragments = ";".join(sorted(BRICS.BRICSDecompose(mol)))
        except Exception:  # noqa: BLE001
            brics_fragments = ""
    rdkit_inchikey = ""
    if compute_rdkit_inchikey:
        try:
            rdkit_inchikey = Chem.MolToInchiKey(mol)
        except Exception:  # noqa: BLE001
            rdkit_inchikey = ""
    return {
        "canonical_smiles": Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        "rdkit_inchikey": rdkit_inchikey,
        "mw": round(float(Descriptors.MolWt(mol)), 4) if compute_physchem_descriptors else "",
        "logp": round(float(Crippen.MolLogP(mol)), 4) if compute_physchem_descriptors else "",
        "tpsa": round(float(rdMolDescriptors.CalcTPSA(mol)), 4) if compute_physchem_descriptors else "",
        "hba": int(Lipinski.NumHAcceptors(mol)) if compute_physchem_descriptors else "",
        "hbd": int(Lipinski.NumHDonors(mol)) if compute_physchem_descriptors else "",
        "rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)) if compute_physchem_descriptors else "",
        "heavy_atoms": int(mol.GetNumHeavyAtoms()),
        "qed": round(float(QED.qed(mol)), 4) if compute_physchem_descriptors else "",
        "murcko_scaffold": scaffold,
        "brics_fragments": brics_fragments,
    }


def build_rank_ready(args: argparse.Namespace, conn: sqlite3.Connection) -> dict[str, Any]:
    path = args.out_dir / "bindingdb_rank_ready_observations.csv"
    counters: Counter[str] = Counter()
    inserted = 0
    descriptor_cache: dict[str, dict[str, Any] | None] = {}
    insert_batch: list[tuple[Any, ...]] = []

    eligible = {
        (row[0], row[1])
        for row in conn.execute(
            """
            SELECT target_id, measurement_type
            FROM aggregated
            GROUP BY target_id, measurement_type
            HAVING COUNT(*) >= ?
            """,
            (args.min_scaffold_group_size,),
        )
    }

    query = """
        SELECT *
        FROM aggregated
        ORDER BY target_id, measurement_type, p_value_median DESC
    """
    conn.row_factory = sqlite3.Row
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RANK_READY_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for row in conn.execute(query):
            if (row["target_id"], row["measurement_type"]) not in eligible:
                continue
            smiles = clean_text(row["smiles_raw"])
            cache_key = row["connectivity_key"] or smiles
            if cache_key not in descriptor_cache:
                descriptor_cache[cache_key] = descriptor_record(
                    smiles,
                    compute_brics_fragments=args.compute_brics_fragments,
                    compute_rdkit_inchikey=args.compute_rdkit_inchikey,
                    compute_physchem_descriptors=args.compute_physchem_descriptors,
                )
            desc = descriptor_cache[cache_key]
            if desc is None or not desc.get("murcko_scaffold"):
                counters["skipped_invalid_or_no_scaffold"] += 1
                continue
            out = {field: row[field] for field in AGGREGATED_FIELDS}
            out.update(desc)
            writer.writerow(out)
            insert_batch.append(tuple(out[field] for field in RANK_READY_FIELDS))
            inserted += 1
            counters[f"{row['measurement_type']}:rank_ready"] += 1
            if len(insert_batch) >= 5000:
                conn.executemany(
                    f"INSERT INTO rank_ready VALUES ({','.join(['?'] * len(RANK_READY_FIELDS))})",
                    insert_batch,
                )
                insert_batch.clear()
            if inserted and inserted % 100000 == 0:
                print(f"[rank_ready] rows={inserted:,}", flush=True)

    if insert_batch:
        conn.executemany(
            f"INSERT INTO rank_ready VALUES ({','.join(['?'] * len(RANK_READY_FIELDS))})",
            insert_batch,
        )
    conn.commit()
    conn.execute("CREATE INDEX IF NOT EXISTS idx_rank_group ON rank_ready(target_id,measurement_type,murcko_scaffold)")
    conn.commit()
    return {"rank_ready_rows": inserted, "counters": dict(sorted(counters.items())), "output": str(path)}


@dataclass
class Candidate:
    aggregate_uid: str
    target_id: str
    target_name: str
    target_organism: str
    target_sequence_hash: str
    target_sequence: str
    measurement_type: str
    measurement_group: str
    smiles: str
    inchikey: str
    connectivity_key: str
    p_value: float
    evidence_count: int
    mw: str
    logp: str
    qed: str
    scaffold: str
    scaffold_heavy_atoms: int
    fp: Any


def mol_fp(smiles: str) -> Any | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)


def heavy_atoms(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    return int(mol.GetNumHeavyAtoms()) if mol is not None else 0


def tanimoto(a: Any, b: Any) -> float:
    if a is None or b is None:
        return 0.0
    return float(DataStructs.TanimotoSimilarity(a, b))


def safe_text_number(value: Any, digits: int = 4) -> str:
    if value in ("", None):
        return ""
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return ""
    if not math.isfinite(numeric):
        return ""
    return str(round(numeric, digits))


def instruction_for(row: Candidate) -> str:
    target = row.target_name or row.target_id
    return (
        f"Given the input molecule and target protein {target}, optimize {row.measurement_type} "
        f"({row.measurement_group}) by proposing a structurally related ligand with stronger measured p-scale activity. "
        "Use the same protein target and measurement type; higher p-value means stronger binding or activity."
    )


def generate_pairs_and_triplets(args: argparse.Namespace, conn: sqlite3.Connection) -> dict[str, Any]:
    triplet_path = args.out_dir / "bindingdb_triplets.csv"
    pair_path = args.out_dir / "bindingdb_rank_pairs.csv"
    all_path = args.out_dir / "all.csv"
    counters: Counter[str] = Counter()
    triplets = 0
    pairs = 0

    target_sequences = {
        row[0]: row[1]
        for row in conn.execute("SELECT target_id,target_sequence FROM targets")
    }
    conn.row_factory = sqlite3.Row
    group_rows = conn.execute(
        """
        SELECT target_id, measurement_type, murcko_scaffold, COUNT(*) AS n
        FROM rank_ready
        GROUP BY target_id, measurement_type, murcko_scaffold
        HAVING n >= ?
        ORDER BY n DESC
        """,
        (args.min_scaffold_group_size,),
    ).fetchall()

    query_count_by_target_metric: Counter[str] = Counter()
    with (
        triplet_path.open("w", encoding="utf-8", newline="") as triplet_handle,
        pair_path.open("w", encoding="utf-8", newline="") as pair_handle,
        all_path.open("w", encoding="utf-8", newline="") as all_handle,
    ):
        triplet_writer = csv.DictWriter(triplet_handle, fieldnames=TRIPLET_FIELDS, extrasaction="ignore")
        pair_writer = csv.DictWriter(pair_handle, fieldnames=PAIR_FIELDS, extrasaction="ignore")
        all_writer = csv.DictWriter(all_handle, fieldnames=TRIPLET_FIELDS, extrasaction="ignore")
        triplet_writer.writeheader()
        pair_writer.writeheader()
        all_writer.writeheader()

        for group_index, group in enumerate(group_rows, start=1):
            target_id = group["target_id"]
            measurement_type = group["measurement_type"]
            scaffold = group["murcko_scaffold"]
            target_metric_key = f"{target_id}|{measurement_type}"
            if args.max_queries_per_target_metric and query_count_by_target_metric[target_metric_key] >= args.max_queries_per_target_metric:
                continue
            records = conn.execute(
                """
                SELECT *
                FROM rank_ready
                WHERE target_id = ? AND measurement_type = ? AND murcko_scaffold = ?
                ORDER BY p_value_median ASC
                """,
                (target_id, measurement_type, scaffold),
            ).fetchall()
            candidates: list[Candidate] = []
            scaffold_atoms = heavy_atoms(scaffold)
            for r in records:
                fp = mol_fp(r["canonical_smiles"])
                if fp is None:
                    continue
                candidates.append(
                    Candidate(
                        aggregate_uid=r["aggregate_uid"],
                        target_id=r["target_id"],
                        target_name=r["target_name"],
                        target_organism=r["target_organism"],
                        target_sequence_hash=r["target_sequence_hash"],
                        target_sequence=target_sequences.get(r["target_id"], ""),
                        measurement_type=r["measurement_type"],
                        measurement_group=r["measurement_group"],
                        smiles=r["canonical_smiles"],
                        inchikey=r["rdkit_inchikey"] or r["inchikey"],
                        connectivity_key=r["connectivity_key"],
                        p_value=float(r["p_value_median"]),
                        evidence_count=int(r["evidence_count"]),
                        mw=safe_text_number(r["mw"], 4),
                        logp=safe_text_number(r["logp"], 4),
                        qed=safe_text_number(r["qed"], 4),
                        scaffold=r["murcko_scaffold"],
                        scaffold_heavy_atoms=scaffold_atoms,
                        fp=fp,
                    )
                )
            if len(candidates) < args.min_scaffold_group_size:
                continue

            for i, input_candidate in enumerate(candidates):
                if args.max_queries_per_target_metric and query_count_by_target_metric[target_metric_key] >= args.max_queries_per_target_metric:
                    break
                positive: Candidate | None = None
                positive_sim = 0.0
                scans = 0
                for cand in reversed(candidates):
                    if cand.connectivity_key == input_candidate.connectivity_key:
                        continue
                    if cand.p_value < input_candidate.p_value + args.min_positive_delta:
                        break
                    scans += 1
                    sim = tanimoto(input_candidate.fp, cand.fp)
                    if sim >= args.min_input_candidate_tanimoto:
                        positive = cand
                        positive_sim = sim
                        break
                    if scans >= args.max_candidate_scan:
                        break
                if positive is None:
                    continue

                negative: Candidate | None = None
                input_negative_sim = 0.0
                positive_negative_sim = 0.0
                scans = 0
                for cand in candidates:
                    if cand.connectivity_key in {input_candidate.connectivity_key, positive.connectivity_key}:
                        continue
                    if cand.p_value > input_candidate.p_value - args.min_negative_delta:
                        break
                    scans += 1
                    sim_in = tanimoto(input_candidate.fp, cand.fp)
                    sim_pos = tanimoto(positive.fp, cand.fp)
                    if sim_in >= args.min_input_candidate_tanimoto and sim_pos >= args.min_positive_negative_tanimoto:
                        negative = cand
                        input_negative_sim = sim_in
                        positive_negative_sim = sim_pos
                        break
                    if scans >= args.max_candidate_scan:
                        break
                if negative is None:
                    continue

                sample_id = f"bindingdbtriplet_{short_hash(input_candidate.aggregate_uid + positive.aggregate_uid + negative.aggregate_uid, 18)}"
                positive_delta = positive.p_value - input_candidate.p_value
                negative_delta = negative.p_value - input_candidate.p_value
                triplet = {
                    "sample_id": sample_id,
                    "source": "BindingDB",
                    "source_version": args.source_version,
                    "task_type": "target_conditioned_binding_ranking_triplet",
                    "instruction": instruction_for(input_candidate),
                    "target_id": input_candidate.target_id,
                    "target_name": input_candidate.target_name,
                    "target_organism": input_candidate.target_organism,
                    "target_sequence_hash": input_candidate.target_sequence_hash,
                    "target_sequence": input_candidate.target_sequence,
                    "measurement_type": input_candidate.measurement_type,
                    "measurement_group": input_candidate.measurement_group,
                    "input_smiles": input_candidate.smiles,
                    "input_inchikey": input_candidate.inchikey,
                    "input_connectivity_key": input_candidate.connectivity_key,
                    "input_p_value": round(input_candidate.p_value, 6),
                    "positive_smiles": positive.smiles,
                    "positive_inchikey": positive.inchikey,
                    "positive_connectivity_key": positive.connectivity_key,
                    "positive_p_value": round(positive.p_value, 6),
                    "positive_delta": round(positive_delta, 6),
                    "negative_smiles": negative.smiles,
                    "negative_inchikey": negative.inchikey,
                    "negative_connectivity_key": negative.connectivity_key,
                    "negative_p_value": round(negative.p_value, 6),
                    "negative_delta": round(negative_delta, 6),
                    "input_positive_tanimoto": round(positive_sim, 6),
                    "input_negative_tanimoto": round(input_negative_sim, 6),
                    "positive_negative_tanimoto": round(positive_negative_sim, 6),
                    "shared_core": scaffold,
                    "shared_core_heavy_atoms": input_candidate.scaffold_heavy_atoms,
                    "input_mw": input_candidate.mw,
                    "positive_mw": positive.mw,
                    "negative_mw": negative.mw,
                    "input_logp": input_candidate.logp,
                    "positive_logp": positive.logp,
                    "negative_logp": negative.logp,
                    "input_qed": input_candidate.qed,
                    "positive_qed": positive.qed,
                    "negative_qed": negative.qed,
                    "input_evidence_count": input_candidate.evidence_count,
                    "positive_evidence_count": positive.evidence_count,
                    "negative_evidence_count": negative.evidence_count,
                }
                pair = {
                    "pair_id": sample_id.replace("bindingdbtriplet_", "bindingdbpair_"),
                    "source": "BindingDB",
                    "source_version": args.source_version,
                    "target_id": input_candidate.target_id,
                    "target_name": input_candidate.target_name,
                    "target_sequence_hash": input_candidate.target_sequence_hash,
                    "measurement_type": input_candidate.measurement_type,
                    "measurement_group": input_candidate.measurement_group,
                    "input_smiles": input_candidate.smiles,
                    "input_connectivity_key": input_candidate.connectivity_key,
                    "input_p_value": round(input_candidate.p_value, 6),
                    "positive_smiles": positive.smiles,
                    "positive_connectivity_key": positive.connectivity_key,
                    "positive_p_value": round(positive.p_value, 6),
                    "positive_delta": round(positive_delta, 6),
                    "input_positive_tanimoto": round(positive_sim, 6),
                    "shared_core": scaffold,
                    "input_evidence_count": input_candidate.evidence_count,
                    "positive_evidence_count": positive.evidence_count,
                }
                triplet_writer.writerow(triplet)
                all_writer.writerow(triplet)
                pair_writer.writerow(pair)
                triplets += 1
                pairs += 1
                query_count_by_target_metric[target_metric_key] += 1
                counters[f"{input_candidate.measurement_type}:triplets"] += 1

            if group_index % 1000 == 0:
                print(f"[triplets] groups={group_index:,}/{len(group_rows):,} triplets={triplets:,}", flush=True)

    return {
        "triplets": triplets,
        "rank_pairs": pairs,
        "counters": dict(sorted(counters.items())),
        "outputs": {
            "all": str(all_path),
            "triplets": str(triplet_path),
            "rank_pairs": str(pair_path),
        },
    }


class DSU:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def find(self, x: str) -> str:
        if x not in self.parent:
            self.parent[x] = x
            return x
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            if stable_hash(ra) > stable_hash(rb):
                ra, rb = rb, ra
            self.parent[rb] = ra


def row_ligands(row: dict[str, str]) -> set[str]:
    return {
        clean_text(row["input_connectivity_key"]),
        clean_text(row["positive_connectivity_key"]),
        clean_text(row["negative_connectivity_key"]),
    }


def split_dataset(args: argparse.Namespace) -> dict[str, Any]:
    all_path = args.out_dir / "all.csv"
    split_paths = {
        "train": args.out_dir / "train.csv",
        "val": args.out_dir / "val.csv",
        "test_seen_target": args.out_dir / "test_seen_target.csv",
        "test_unseen_target": args.out_dir / "test_unseen_target.csv",
    }
    with all_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    target_ids = sorted({row["target_id"] for row in rows})
    unseen_targets = {target for target in target_ids if stable_hash(target) % 100 < args.unseen_target_bucket}
    unseen_rows = [row for row in rows if row["target_id"] in unseen_targets]
    unseen_ligands: set[str] = set()
    for row in unseen_rows:
        unseen_ligands.update(row_ligands(row))

    seen_rows = [
        row
        for row in rows
        if row["target_id"] not in unseen_targets and row_ligands(row).isdisjoint(unseen_ligands)
    ]
    dsu = DSU()
    for row in seen_rows:
        ligands = list(row_ligands(row))
        if not ligands:
            continue
        for ligand in ligands[1:]:
            dsu.union(ligands[0], ligand)

    splits: dict[str, list[dict[str, str]]] = {
        "train": [],
        "val": [],
        "test_seen_target": [],
        "test_unseen_target": unseen_rows,
    }
    for row in seen_rows:
        root = dsu.find(row["input_connectivity_key"])
        bucket = stable_hash(root) % 100
        if bucket < 100 - args.seen_val_bucket - args.seen_test_bucket:
            splits["train"].append(row)
        elif bucket < 100 - args.seen_test_bucket:
            splits["val"].append(row)
        else:
            splits["test_seen_target"].append(row)

    for name, path in split_paths.items():
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=TRIPLET_FIELDS, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(splits[name])

    def ligand_set(name: str) -> set[str]:
        out: set[str] = set()
        for row in splits[name]:
            out.update(row_ligands(row))
        return out

    def target_set(name: str) -> set[str]:
        return {row["target_id"] for row in splits[name]}

    leakage = {}
    split_names = list(splits)
    for i, a in enumerate(split_names):
        for b in split_names[i + 1 :]:
            leakage[f"ligand_overlap_{a}_{b}"] = len(ligand_set(a) & ligand_set(b))
    for name in ["train", "val", "test_seen_target"]:
        leakage[f"target_overlap_{name}_test_unseen_target"] = len(target_set(name) & target_set("test_unseen_target"))

    return {
        "rows": {name: len(value) for name, value in splits.items()},
        "unique_ligands": {name: len(ligand_set(name)) for name in splits},
        "unique_targets": {name: len(target_set(name)) for name in splits},
        "dropped_seen_rows_due_to_unseen_ligand_overlap": len(rows) - len(unseen_rows) - len(seen_rows),
        "leakage": leakage,
        "outputs": {name: str(path) for name, path in split_paths.items()},
    }


def write_metric_summary(args: argparse.Namespace, conn: sqlite3.Connection) -> dict[str, Any]:
    path = args.out_dir / "bindingdb_metric_summary.csv"
    fields = [
        "measurement_group",
        "measurement_type",
        "aggregated_ligand_target_rows",
        "unique_ligands",
        "unique_targets",
        "median_p_value",
    ]
    rows = []
    for group, mtype, n, ligands, targets in conn.execute(
        """
        SELECT measurement_group, measurement_type, COUNT(*), COUNT(DISTINCT connectivity_key), COUNT(DISTINCT target_id)
        FROM aggregated
        GROUP BY measurement_group, measurement_type
        ORDER BY measurement_group, measurement_type
        """
    ):
        values = [
            r[0]
            for r in conn.execute(
                "SELECT p_value_median FROM aggregated WHERE measurement_group=? AND measurement_type=?",
                (group, mtype),
            )
        ]
        rows.append([group, mtype, n, ligands, targets, round(float(median(values)), 6) if values else ""])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        writer.writerows(rows)
    return {"rows": len(rows), "output": str(path)}


def write_readme(args: argparse.Namespace, stats: dict[str, Any]) -> None:
    path = args.out_dir / "README.md"
    text = f"""# BindingDB Target-Conditioned Ranking Dataset

This folder was generated from the BindingDB full TSV snapshot.

## Task

`input ligand SMILES + target protein context + measurement instruction -> stronger positive ligand + weaker hard-negative ligand`

The final supervised rows are target-conditioned triplets. Each row compares ligands only within the same target protein and the same measurement type.

## Measurement normalization

- Ki and Kd are assigned to `binding_affinity`.
- IC50 and EC50 are assigned to `activity_potency`.
- All nM values are converted to p-scale by `p = 9 - log10(value_nM)`.
- Higher p-scale values mean stronger binding/activity.
- Duplicate ligand-target-measurement records are aggregated by median p-scale.
- Censored relations such as `<`, `>`, `<=`, `>=` are kept in the raw observation file but excluded from ranking/triplet construction.
- BindingDB full TSV does not expose a reliable agonist/antagonist column, so `modulation_numeric=0.5` means unknown, not neutral.

## Output files

- `bindingdb_observations.csv`: parsed Ki/Kd/IC50/EC50 observations.
- `bindingdb_aggregated_observations.csv`: median p-scale per ligand + target + measurement type.
- `bindingdb_rank_ready_observations.csv`: aggregate rows with RDKit canonical SMILES and scaffold. MW/LogP/QED descriptor columns are filled only when `--compute-physchem-descriptors` is enabled.
- `bindingdb_rank_pairs.csv`: input-positive ranking pairs.
- `bindingdb_triplets.csv` and `all.csv`: input-positive-negative target-conditioned triplets.
- `train.csv`, `val.csv`, `test_seen_target.csv`, `test_unseen_target.csv`: leakage-controlled splits.
- `bindingdb_target_summary.csv`: target/protein metadata and sequences.
- `bindingdb_metric_summary.csv`: counts by measurement group/type.
- `dataset_stats.json`: full processing statistics.

## Current stats

```json
{json.dumps(stats, indent=2, ensure_ascii=False)}
```
"""
    path.write_text(text, encoding="utf-8")


def main() -> int:
    safe_field_size_limit()
    args = parse_args()
    if not args.input.is_file():
        raise SystemExit(f"ERROR: BindingDB input not found: {args.input}")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    db_path = args.out_dir / "bindingdb_pipeline.sqlite"
    conn = init_db(db_path)

    started = datetime.now(timezone.utc)
    stats: dict[str, Any] = {
        "created_at_utc": started.isoformat(),
        "source": {
            "name": "BindingDB",
            "source_version": args.source_version,
            "input": str(args.input),
        },
        "parameters": {
            "max_raw_rows": args.max_raw_rows,
            "min_positive_delta": args.min_positive_delta,
            "min_negative_delta": args.min_negative_delta,
            "min_input_candidate_tanimoto": args.min_input_candidate_tanimoto,
            "min_positive_negative_tanimoto": args.min_positive_negative_tanimoto,
            "min_scaffold_group_size": args.min_scaffold_group_size,
            "max_candidate_scan": args.max_candidate_scan,
            "max_queries_per_target_metric": args.max_queries_per_target_metric,
            "compute_brics_fragments": args.compute_brics_fragments,
            "compute_rdkit_inchikey": args.compute_rdkit_inchikey,
            "compute_physchem_descriptors": args.compute_physchem_descriptors,
            "split": {
                "seen_val_bucket": args.seen_val_bucket,
                "seen_test_bucket": args.seen_test_bucket,
                "unseen_target_bucket": args.unseen_target_bucket,
            },
        },
    }
    try:
        print("[pipeline] normalize", flush=True)
        stats["normalization"] = normalize_to_observations(args, conn)
        print("[pipeline] target summary", flush=True)
        stats["target_summary"] = write_target_summary(args, conn)
        print("[pipeline] aggregate", flush=True)
        stats["aggregation"] = aggregate_exact_observations(args, conn)
        print("[pipeline] metric summary", flush=True)
        stats["metric_summary"] = write_metric_summary(args, conn)
        print("[pipeline] rank-ready descriptors", flush=True)
        stats["rank_ready"] = build_rank_ready(args, conn)
        print("[pipeline] pairs/triplets", flush=True)
        stats["ranking"] = generate_pairs_and_triplets(args, conn)
        print("[pipeline] split", flush=True)
        stats["split"] = split_dataset(args)
        stats["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
        stats["duration_seconds"] = round((datetime.now(timezone.utc) - started).total_seconds(), 2)
        stats_path = args.out_dir / "dataset_stats.json"
        stats_path.write_text(json.dumps(stats, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        write_readme(args, stats)
        print(json.dumps({"dataset_stats": str(stats_path), "split_rows": stats["split"]["rows"]}, indent=2), flush=True)
    finally:
        conn.close()
        if not args.keep_sqlite and db_path.exists():
            db_path.unlink()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
