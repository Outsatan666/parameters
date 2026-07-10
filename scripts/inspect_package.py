from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from scripts.mol2 import Mol2ParseError, parse_mol2, sha256_file

ISSUE_TOKENS = ("warning", "error", "failed", "missing", "unknown", "unsupported")
CHARGE_TOLERANCE = 1e-4
PARAMETER_SUFFIXES = (".str", ".par", ".prm", ".rtf")


@dataclass(frozen=True)
class Issue:
    code: str
    severity: str
    message: str
    context: str | None = None


def inventory_files(root: str | Path) -> list[str]:
    root_path = Path(root)
    return sorted(str(path.relative_to(root_path)).replace("\\", "/") for path in root_path.rglob("*") if path.is_file())


def classify_branch(path: str | Path) -> str | None:
    normalized = str(path).replace("\\", "/").lower()
    name = Path(path).name.lower()
    if "match" in normalized or name.startswith("match") or ".match." in name:
        return "MATCH"
    if "mmff" in normalized or name.startswith("mmff") or ".mmff." in name:
        return "MMFF"
    return None


def branch_inventory(root: str | Path) -> dict[str, list[str]]:
    result = {"MATCH": [], "MMFF": [], "UNCLASSIFIED": []}
    for relative in inventory_files(root):
        branch = classify_branch(relative)
        result[branch or "UNCLASSIFIED"].append(relative)
    return result


def scan_logs(root: str | Path, context_radius: int = 1) -> list[Issue]:
    root_path = Path(root)
    issues: list[Issue] = []
    for path in root_path.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in {".log", ".out", ".err", ".txt"}:
            continue
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        for index, line in enumerate(lines):
            lowered = line.lower()
            matched = next((token for token in ISSUE_TOKENS if token in lowered), None)
            if matched is None:
                continue
            start = max(0, index - context_radius)
            end = min(len(lines), index + context_radius + 1)
            severity = "CRITICAL" if matched in {"error", "failed"} and "0 error" not in lowered else "REVIEW"
            issues.append(Issue(f"LOG_{matched.upper()}", severity, f"{path.name}:{index + 1}: {line.strip()}", "\n".join(lines[start:end])))
    return issues


def _output_mol2_candidates(root: Path, branch: str) -> list[Path]:
    return sorted(path for path in root.rglob("*.mol2") if classify_branch(path.relative_to(root)) == branch)


def run_qc(input_mol2: str | Path, extracted_root: str | Path, *, expected_charge: int | None = None) -> dict:
    input_path = Path(input_mol2)
    root = Path(extracted_root)
    issues: list[Issue] = []
    try:
        input_data = parse_mol2(input_path)
    except Mol2ParseError as exc:
        return {"status": "FAILED", "input_sha256": sha256_file(input_path), "issues": [asdict(Issue("INPUT_MALFORMED", "CRITICAL", str(exc)))]}

    if input_data.duplicate_atom_names:
        issues.append(Issue("DUPLICATE_ATOM_NAMES", "REVIEW", "Input MOL2 has non-unique atom names; reliable name-based output mapping may be impossible: " + ", ".join(input_data.duplicate_atom_names)))

    branches = branch_inventory(root)
    if not branches["MATCH"]:
        issues.append(Issue("MATCH_MISSING", "CRITICAL", "MATCH branch was not identified"))
    elif not any(Path(name).suffix.lower() in PARAMETER_SUFFIXES for name in branches["MATCH"]):
        issues.append(Issue("MATCH_PARAMETERS_MISSING", "CRITICAL", "MATCH branch has no CHARMM parameter files (.str/.par/.prm/.rtf); found only: " + ", ".join(branches["MATCH"])))

    charge_checks: dict[str, dict] = {}
    atom_counts: dict[str, int | None] = {"INPUT": len(input_data.atoms)}
    for branch in ("MATCH", "MMFF"):
        candidates = _output_mol2_candidates(root, branch)
        if not candidates:
            atom_counts[branch] = None
            continue
        try:
            output_data = parse_mol2(candidates[0])
        except Mol2ParseError as exc:
            issues.append(Issue(f"{branch}_MOL2_PARSE_FAILED", "CRITICAL", f"{candidates[0].name}: {exc}"))
            atom_counts[branch] = None
            continue
        atom_counts[branch] = len(output_data.atoms)
        if len(output_data.atoms) != len(input_data.atoms):
            issues.append(Issue(f"{branch}_ATOM_COUNT_MISMATCH", "REVIEW", f"Input atom count {len(input_data.atoms)} != {branch} atom count {len(output_data.atoms)}"))
        if expected_charge is not None and output_data.total_charge is not None:
            delta = output_data.total_charge - expected_charge
            charge_checks[branch] = {"sum": output_data.total_charge, "expected": expected_charge, "delta": delta}
            if abs(delta) > CHARGE_TOLERANCE:
                issues.append(Issue(f"{branch}_CHARGE_MISMATCH", "REVIEW", f"{branch} charge sum {output_data.total_charge:.8f} differs from expected {expected_charge} by {delta:.8f}"))

    issues.extend(scan_logs(root))
    status = "REVIEW" if issues else "TECHNICAL_PASS"
    return {"status": status, "input_sha256": sha256_file(input_path), "input": {"name": input_data.name, "atom_count": len(input_data.atoms), "bond_count": len(input_data.bonds), "charge_type": input_data.charge_type, "input_charge_sum": input_data.total_charge, "duplicate_atom_names": list(input_data.duplicate_atom_names)}, "branches": branches, "atom_counts": atom_counts, "charge_checks": charge_checks, "issues": [asdict(issue) for issue in issues]}


def write_qc(qc: dict, qc_json: str | Path, issues_tsv: str | Path) -> None:
    qc_path = Path(qc_json)
    qc_path.parent.mkdir(parents=True, exist_ok=True)
    qc_path.write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    lines = ["code\tseverity\tmessage\tcontext"]
    for issue in qc.get("issues", []):
        values = [str(issue.get("code", "")), str(issue.get("severity", "")), str(issue.get("message", "")).replace("\t", " ").replace("\n", "\\n"), str(issue.get("context", "") or "").replace("\t", " ").replace("\n", "\\n")]
        lines.append("\t".join(values))
    Path(issues_tsv).write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Technical QC for a SwissParam package")
    parser.add_argument("input_mol2", type=Path)
    parser.add_argument("extracted_root", type=Path)
    parser.add_argument("--expected-charge", type=int)
    parser.add_argument("--qc-json", type=Path, default=Path("qc.json"))
    parser.add_argument("--issues-tsv", type=Path, default=Path("issues.tsv"))
    args = parser.parse_args()
    qc = run_qc(args.input_mol2, args.extracted_root, expected_charge=args.expected_charge)
    write_qc(qc, args.qc_json, args.issues_tsv)
    print(qc["status"])
    return 0 if qc["status"] != "FAILED" else 1


if __name__ == "__main__":
    raise SystemExit(main())
