from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


def build_package(molecule: str, *, original_mol2: str | Path, raw_archive: str | Path, extracted_root: str | Path, comparison_tsv: str | Path, comparison_md: str | Path, qc_json: str | Path, issues_tsv: str | Path, manifest_entry: dict, packages_dir: str | Path = "packages", staging_root: str | Path = ".package_staging") -> Path:
    staging_base = Path(staging_root)
    molecule_root = staging_base / molecule
    if molecule_root.exists():
        shutil.rmtree(molecule_root)
    for rel in ("input", "swissparam/raw", "swissparam/MATCH", "swissparam/MMFF", "comparison", "qc"):
        (molecule_root / rel).mkdir(parents=True, exist_ok=True)
    shutil.copy2(original_mol2, molecule_root / "input/original.mol2")
    shutil.copy2(raw_archive, molecule_root / "swissparam/raw/results.tar.gz")
    shutil.copy2(comparison_tsv, molecule_root / "comparison/comparison.tsv")
    shutil.copy2(comparison_md, molecule_root / "comparison/comparison.md")
    shutil.copy2(qc_json, molecule_root / "qc/qc.json")
    shutil.copy2(issues_tsv, molecule_root / "qc/issues.tsv")
    (molecule_root / "manifest.json").write_text(json.dumps(manifest_entry, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    extracted = Path(extracted_root)
    for path in extracted.rglob("*"):
        if not path.is_file():
            continue
        normalized = str(path.relative_to(extracted)).replace("\\", "/").lower()
        branch = "MATCH" if "match" in normalized else "MMFF" if "mmff" in normalized else None
        if branch:
            target = molecule_root / "swissparam" / branch / path.name
            if target.exists():
                target = target.with_name(f"{path.parent.name}_{path.name}")
            shutil.copy2(path, target)
    output_dir = Path(packages_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / f"{molecule}_CHARMM_STARTING_PARAMETERS.zip"
    with ZipFile(output, "w", compression=ZIP_DEFLATED) as archive:
        for path in sorted(molecule_root.rglob("*")):
            if path.is_file():
                archive.write(path, arcname=str(path.relative_to(staging_base)))
    return output


def main() -> int:
    p = argparse.ArgumentParser(description="Build per-molecule parameter ZIP")
    p.add_argument("molecule")
    for name in ("original-mol2", "raw-archive", "extracted-root", "comparison-tsv", "comparison-md", "qc-json", "issues-tsv", "manifest-entry-json"):
        p.add_argument(f"--{name}", required=True)
    a = p.parse_args()
    manifest_entry = json.loads(Path(a.manifest_entry_json).read_text(encoding="utf-8"))
    print(build_package(a.molecule, original_mol2=a.original_mol2, raw_archive=a.raw_archive, extracted_root=a.extracted_root, comparison_tsv=a.comparison_tsv, comparison_md=a.comparison_md, qc_json=a.qc_json, issues_tsv=a.issues_tsv, manifest_entry=manifest_entry))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
