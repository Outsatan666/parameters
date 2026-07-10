from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

from scripts.mol2 import parse_mol2
from scripts.qm_inputs import build_crest_command, build_crest_xyz, build_orca_input
from scripts.structure_files import StructureFileError, build_pdb, generate_psf, psfgen_available

DEFAULT_ORCA_KEYWORDS = "! R2SCAN-3c OPT freq"


def write_derived_artifacts(
    mol2_path: str | Path,
    out_dir: str | Path,
    *,
    orca_keywords: str = DEFAULT_ORCA_KEYWORDS,
    charge: int = 0,
    multiplicity: int = 1,
    solvent: str | None = None,
    topology_paths: list[str] | None = None,
    segid: str = "MOL",
) -> dict[str, Path]:
    """Generate the FFtk-consistent PDB and the ORCA/CREST inputs for one molecule.

    A PSF is produced only when CHARMM topology is supplied *and* psfgen/VMD is on
    PATH; otherwise it is skipped (never faked). Returns a map of artifact name -> path.
    """
    mol2_path = Path(mol2_path)
    mol2 = parse_mol2(mol2_path)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    stem = mol2_path.stem
    written: dict[str, Path] = {}

    pdb_path = out / f"{stem}.pdb"
    pdb_path.write_text(build_pdb(mol2), encoding="utf-8")
    written["pdb"] = pdb_path

    orca_path = out / f"{stem}_orca.inp"
    orca_path.write_text(build_orca_input(mol2, keywords=orca_keywords, charge=charge, multiplicity=multiplicity, solvent=solvent), encoding="utf-8")
    written["orca"] = orca_path

    xyz_path = out / f"{stem}.xyz"
    xyz_path.write_text(build_crest_xyz(mol2), encoding="utf-8")
    written["crest_xyz"] = xyz_path

    crest_cmd_path = out / f"{stem}_crest.sh"
    crest_cmd_path.write_text(build_crest_command(xyz_path.name, solvent=solvent.lower() if solvent else None, charge=charge) + "\n", encoding="utf-8")
    written["crest_cmd"] = crest_cmd_path

    if topology_paths and psfgen_available():
        try:
            written["psf"] = generate_psf(topology_paths, str(pdb_path), out, segid=segid, output_prefix=stem)
        except (StructureFileError, OSError, subprocess.SubprocessError):
            # psfgen present but failed (e.g. bad/absent topology): keep PDB+QM, skip PSF.
            pass

    return written


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate PDB + ORCA/CREST inputs for one MOL2")
    parser.add_argument("mol2", type=Path)
    parser.add_argument("--out-dir", type=Path, default=Path("derived"))
    parser.add_argument("--orca-keywords", default=DEFAULT_ORCA_KEYWORDS)
    parser.add_argument("--charge", type=int, default=0)
    parser.add_argument("--multiplicity", type=int, default=1)
    parser.add_argument("--solvent")
    parser.add_argument("--topology", action="append", default=[])
    parser.add_argument("--segid", default="MOL")
    args = parser.parse_args()
    written = write_derived_artifacts(args.mol2, args.out_dir, orca_keywords=args.orca_keywords, charge=args.charge, multiplicity=args.multiplicity, solvent=args.solvent, topology_paths=args.topology or None, segid=args.segid)
    for label, path in written.items():
        print(f"{label}\t{path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
