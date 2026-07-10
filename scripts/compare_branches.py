from __future__ import annotations

import argparse
import csv
from pathlib import Path

from scripts.inspect_package import classify_branch
from scripts.mol2 import parse_mol2


def _find_branch_mol2(root: Path, branch: str) -> Path | None:
    candidates = [p for p in root.rglob("*.mol2") if classify_branch(p.relative_to(root)) == branch]
    return sorted(candidates)[0] if candidates else None


def compare_branches(extracted_root: str | Path) -> dict:
    root = Path(extracted_root)
    match_path = _find_branch_mol2(root, "MATCH")
    mmff_path = _find_branch_mol2(root, "MMFF")
    result = {"match_path": str(match_path) if match_path else None, "mmff_path": str(mmff_path) if mmff_path else None, "mapping_mode": None, "rows": [], "summary": {}}
    if match_path is None or mmff_path is None:
        result["summary"] = {"comparable": False, "reason": "MATCH or MMFF MOL2 missing"}
        return result
    match = parse_mol2(match_path)
    mmff = parse_mol2(mmff_path)
    match_names = [a.name for a in match.atoms]
    mmff_names = [a.name for a in mmff.atoms]
    unique_names = len(match_names) == len(set(match_names)) and len(mmff_names) == len(set(mmff_names)) and set(match_names) == set(mmff_names)
    if unique_names:
        result["mapping_mode"] = "atom_name"
        mmff_by_name = {a.name: a for a in mmff.atoms}
        pairs = [(a, mmff_by_name[a.name]) for a in match.atoms]
    elif len(match.atoms) == len(mmff.atoms):
        result["mapping_mode"] = "index_fallback_unverified"
        pairs = list(zip(match.atoms, mmff.atoms, strict=True))
    else:
        result["summary"] = {"comparable": False, "reason": "Atom counts differ and unique name mapping is unavailable"}
        return result
    rows = []
    for ma, mm in pairs:
        delta = None if ma.charge is None or mm.charge is None else ma.charge - mm.charge
        rows.append({"match_index": ma.atom_id, "match_name": ma.name, "match_type": ma.atom_type, "match_charge": ma.charge, "mmff_index": mm.atom_id, "mmff_name": mm.name, "mmff_type": mm.atom_type, "mmff_charge": mm.charge, "charge_delta_match_minus_mmff": delta, "atom_type_equal": ma.atom_type == mm.atom_type})
    result["rows"] = rows
    result["summary"] = {"comparable": True, "atom_count_match": len(match.atoms), "atom_count_mmff": len(mmff.atoms), "match_total_charge": match.total_charge, "mmff_total_charge": mmff.total_charge, "different_atom_types": sum(not row["atom_type_equal"] for row in rows), "max_abs_charge_delta": max((abs(row["charge_delta_match_minus_mmff"]) for row in rows if row["charge_delta_match_minus_mmff"] is not None), default=None)}
    return result


def write_comparison(result: dict, tsv_path: str | Path, md_path: str | Path) -> None:
    tsv, md = Path(tsv_path), Path(md_path)
    tsv.parent.mkdir(parents=True, exist_ok=True)
    md.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["match_index", "match_name", "match_type", "match_charge", "mmff_index", "mmff_name", "mmff_type", "mmff_charge", "charge_delta_match_minus_mmff", "atom_type_equal"]
    with tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(result.get("rows", []))
    summary = result.get("summary", {})
    lines = ["# MATCH vs MMFF comparison", "", f"- Mapping mode: `{result.get('mapping_mode')}`", f"- Comparable: `{summary.get('comparable')}`"]
    lines.extend(f"- {k}: `{v}`" for k, v in summary.items() if k != "comparable")
    lines.extend(["", "> MATCH and MMFF are independent starting branches. No parameter averaging or automatic branch mixing is performed.", ""])
    md.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    p = argparse.ArgumentParser(description="Compare SwissParam MATCH and MMFF branches")
    p.add_argument("extracted_root", type=Path)
    p.add_argument("--tsv", type=Path, required=True)
    p.add_argument("--md", type=Path, required=True)
    a = p.parse_args()
    write_comparison(compare_branches(a.extracted_root), a.tsv, a.md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
