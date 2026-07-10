from __future__ import annotations

import argparse
import csv
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

from scripts.build_archives import build_package
from scripts.compare_branches import compare_branches, write_comparison
from scripts.inspect_package import run_qc, write_qc
from scripts.manifest import ManifestStore, make_job_key
from scripts.mol2 import Mol2ParseError, parse_mol2, sha256_file
from scripts.swissparam_client import SwissParamClient, SwissParamError, safe_extract_tar


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def load_states(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open(encoding="utf-8", newline="") as handle:
        return {row["molecule"]: row for row in csv.DictReader(handle, delimiter="\t") if row.get("molecule")}


def _expected_charge(row: dict[str, str] | None) -> int | None:
    if not row:
        return None
    value = (row.get("expected_charge") or "").strip()
    return int(value) if value else None


def process_one(mol2_path: Path, *, client: SwissParamClient, manifest: ManifestStore, state_row: dict[str, str] | None, approach: str, poll_interval_s: float, max_total_wait_s: float, force: bool, workspace: Path) -> dict:
    molecule = mol2_path.stem
    input_sha = sha256_file(mol2_path)
    job_key = make_job_key(input_sha, approach=approach)
    existing = manifest.get(job_key)
    if existing and manifest.is_complete(job_key) and not force:
        return {"molecule": molecule, "action": "SKIP_ALREADY_COMPLETE", **existing}
    try:
        input_data = parse_mol2(mol2_path)
    except Mol2ParseError as exc:
        entry = {"molecule": molecule, "status": "FAILED", "review_reason": f"INPUT_MALFORMED: {exc}", "input_sha256": input_sha, "updated_utc": utc_now()}
        manifest.upsert(job_key, **entry)
        return {"molecule": molecule, "action": "FAILED_INPUT", **entry}
    raw_dir = workspace / "results/raw" / molecule
    extracted_dir = raw_dir / "extracted"
    raw_dir.mkdir(parents=True, exist_ok=True)
    session_number = existing.get("session_number") if existing else None
    if session_number and not force:
        action = "RESUME_POLLING"
    else:
        client.check_health()
        session_number = client.submit(mol2_path, approach=approach, raw_response_path=raw_dir / "submit_response.txt")
        action = "SUBMITTED"
        manifest.upsert(job_key, molecule=molecule, status="GENERATED", session_number=session_number, input_sha256=input_sha, generator="swissparam", approach=approach, updated_utc=utc_now(), expected_charge=_expected_charge(state_row) if state_row else "EXPECTED_CHARGE_UNKNOWN")
    status_result = client.poll(session_number, poll_interval_s=poll_interval_s, max_total_wait_s=max_total_wait_s, raw_status_dir=raw_dir / "status")
    if status_result.state == "failed":
        entry = {"molecule": molecule, "status": "FAILED", "session_number": session_number, "review_reason": "SwissParam reported failed", "updated_utc": utc_now()}
        manifest.upsert(job_key, **entry)
        return {"molecule": molecule, "action": action, **entry}
    archive = client.retrieve(session_number, raw_dir / "results.tar.gz")
    if extracted_dir.exists():
        shutil.rmtree(extracted_dir)
    safe_extract_tar(archive, extracted_dir)
    expected_charge = _expected_charge(state_row)
    qc = run_qc(mol2_path, extracted_dir, expected_charge=expected_charge)
    qc_json, issues_tsv = raw_dir / "qc.json", raw_dir / "issues.tsv"
    write_qc(qc, qc_json, issues_tsv)
    comparison = compare_branches(extracted_dir)
    comparison_tsv = workspace / "reports" / f"{molecule}_comparison.tsv"
    comparison_md = workspace / "reports" / f"{molecule}_comparison.md"
    write_comparison(comparison, comparison_tsv, comparison_md)
    entry = {"molecule": molecule, "status": qc["status"], "session_number": session_number, "input_sha256": input_sha, "generator": "swissparam", "approach": approach, "updated_utc": utc_now(), "expected_charge": expected_charge if expected_charge is not None else "EXPECTED_CHARGE_UNKNOWN", "input_atom_count": len(input_data.atoms), "input_bond_count": len(input_data.bonds), "archive_sha256": sha256_file(archive), "issues": qc.get("issues", [])}
    manifest.upsert(job_key, **entry)
    package = build_package(molecule, original_mol2=mol2_path, raw_archive=archive, extracted_root=extracted_dir, comparison_tsv=comparison_tsv, comparison_md=comparison_md, qc_json=qc_json, issues_tsv=issues_tsv, manifest_entry=entry, packages_dir=workspace / "packages", staging_root=workspace / ".package_staging")
    return {"molecule": molecule, "action": action, "package": str(package), **entry}


def write_review_queue(results: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = ["molecule\tstate\tstatus\treview_reason\tmatch_available\tmmff_available\ttier2_recommended\ttier2_backend"]
    for result in results:
        if result.get("status") == "REVIEW":
            issue_codes = ",".join(issue.get("code", "") for issue in result.get("issues", []))
            rows.append("\t".join([str(result.get("molecule", "")), str(result.get("state", "")), "REVIEW", issue_codes, "unknown", "unknown", "yes", "FFParam discovery required"]))
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def main() -> int:
    p = argparse.ArgumentParser(description="Batch SwissParam parameterization")
    p.add_argument("--input-dir", type=Path, default=Path("input"))
    p.add_argument("--states", type=Path, default=Path("config/states.tsv"))
    p.add_argument("--manifest", type=Path, default=Path("manifests/manifest.json"))
    p.add_argument("--approach", choices=["both", "match", "mmff-based"], default="both")
    p.add_argument("--poll-interval", type=float, default=15.0)
    p.add_argument("--max-total-wait", type=float, default=7200.0)
    p.add_argument("--force", action="store_true")
    a = p.parse_args()
    workspace = Path.cwd()
    inputs = sorted(a.input_dir.glob("*.mol2"))
    if not inputs:
        print(f"No MOL2 inputs found in {a.input_dir}")
        return 2
    states = load_states(a.states)
    manifest = ManifestStore.load(a.manifest)
    client = SwissParamClient()
    results = []
    for mol2_path in inputs:
        try:
            result = process_one(mol2_path, client=client, manifest=manifest, state_row=states.get(mol2_path.stem), approach=a.approach, poll_interval_s=a.poll_interval, max_total_wait_s=a.max_total_wait, force=a.force, workspace=workspace)
        except (SwissParamError, OSError, ValueError) as exc:
            input_sha = sha256_file(mol2_path)
            job_key = make_job_key(input_sha, approach=a.approach)
            manifest.upsert(job_key, molecule=mol2_path.stem, status="FAILED", review_reason=str(exc), updated_utc=utc_now())
            result = {"molecule": mol2_path.stem, "status": "FAILED", "action": "EXCEPTION", "review_reason": str(exc)}
        results.append(result)
        print(json.dumps(result, sort_keys=True))
    write_review_queue(results, workspace / "reports/review_queue.tsv")
    counts = {status: sum(r.get("status") == status for r in results) for status in ("GENERATED", "TECHNICAL_PASS", "REVIEW", "FAILED", "SECOND_LEVEL_COMPLETE")}
    print(json.dumps({"total_inputs": len(results), "generated": counts["GENERATED"], "technical_pass": counts["TECHNICAL_PASS"], "review": counts["REVIEW"], "failed": counts["FAILED"], "tier2_complete": counts["SECOND_LEVEL_COMPLETE"]}, sort_keys=True))
    return 1 if counts["FAILED"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
