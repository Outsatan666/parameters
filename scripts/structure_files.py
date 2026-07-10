from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path

from scripts.mol2 import Mol2Data, parse_mol2
from scripts.qm_inputs import element_symbol

# psfgen ships inside VMD; accept either a standalone psfgen or VMD as the driver.
PSFGEN_CANDIDATES = ("psfgen", "vmd")


class StructureFileError(RuntimeError):
    """Raised for PDB/PSF generation problems (bad inputs or missing psfgen)."""


def build_pdb(mol2: Mol2Data, *, residue_name: str | None = None, hetatm: bool = True) -> str:
    """Emit a PDB whose atom serials and names mirror the input MOL2 exactly.

    The atom serial equals the MOL2 ``atom_id`` and the atom name is preserved
    verbatim, so PDB/PSF/MOL2 share one numbering that VMD FFtk can read.
    """
    if not mol2.atoms:
        raise StructureFileError("Cannot build PDB for a molecule with no atoms")
    record = "HETATM" if hetatm else "ATOM  "
    lines: list[str] = []
    for atom in mol2.atoms:
        if len(atom.name) > 4:
            raise StructureFileError(f"Atom name {atom.name!r} exceeds PDB's 4-character field; cannot preserve numbering. Rename the atom or route to REVIEW.")
        resname = (residue_name or atom.subst_name or "MOL")[:3]
        element = element_symbol(atom)
        name = atom.name if len(atom.name) >= 4 else f" {atom.name}"
        lines.append(
            f"{record}{atom.atom_id:>5d} {name:<4.4s}"
            f" {resname:>3s} A{1:>4d}    "
            f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
            f"{1.0:6.2f}{0.0:6.2f}          {element:>2s}"
        )
    lines.append("END")
    return "\n".join(lines) + "\n"


def numbering_matches(mol2: Mol2Data, pdb_text: str) -> bool:
    """True iff the PDB's atom serial/name sequence equals the MOL2's, in order."""
    atom_lines = [ln for ln in pdb_text.splitlines() if ln.startswith(("ATOM", "HETATM"))]
    if len(atom_lines) != len(mol2.atoms):
        return False
    for atom, line in zip(mol2.atoms, atom_lines, strict=True):
        try:
            serial = int(line[6:11])
        except ValueError:
            return False
        if serial != atom.atom_id or line[12:16].strip() != atom.name:
            return False
    return True


def build_psfgen_script(topology_paths: list[str], pdb_path: str, *, segid: str = "MOL", output_prefix: str = "molecule") -> str:
    """Build a psfgen ``.pgn`` script that assembles the PSF from topology + PDB.

    This is the canonical VMD/FFtk route: the same PDB is fed back so the written
    PSF and PDB share identical atom ordering.
    """
    if not topology_paths:
        raise StructureFileError("At least one CHARMM topology (.rtf/.str) is required for psfgen")
    lines = [f"topology {path}" for path in topology_paths]
    lines += [
        f"segment {segid} {{",
        f"    pdb {pdb_path}",
        "    first none",
        "    last none",
        "}",
        f"coordpdb {pdb_path} {segid}",
        "guesscoord",
        f"writepsf {output_prefix}.psf",
        f"writepdb {output_prefix}.pdb",
    ]
    return "\n".join(lines) + "\n"


def psfgen_available() -> str | None:
    """Return the path to psfgen (or VMD, which bundles it), or None if absent.

    Mirrors the repository's inventory pattern: probe only, never install.
    """
    for candidate in PSFGEN_CANDIDATES:
        found = shutil.which(candidate)
        if found:
            return found
    return None


def generate_psf(topology_paths: list[str], pdb_path: str, workdir: str | Path, *, segid: str = "MOL", output_prefix: str = "molecule") -> Path:
    """Run psfgen to produce a PSF. Raises if psfgen/VMD is not installed.

    We do not fabricate a PSF in pure Python: a correct CHARMM PSF needs masses and
    generated angles/dihedrals/impropers from the topology, which is exactly psfgen's
    job and must be validated in VMD FFtk.
    """
    driver = psfgen_available()
    if driver is None:
        raise StructureFileError("psfgen/VMD not found on PATH; PSF generation requires VMD (bundles psfgen). Install VMD or run on a workstation that has it.")
    work = Path(workdir)
    work.mkdir(parents=True, exist_ok=True)
    script_path = work / f"{output_prefix}.pgn"
    script_path.write_text(build_psfgen_script(topology_paths, pdb_path, segid=segid, output_prefix=str(work / output_prefix)), encoding="utf-8")
    command = [driver, "-dispdev", "text", "-e", str(script_path)] if Path(driver).stem.lower() == "vmd" else [driver, str(script_path)]
    try:
        subprocess.run(command, check=True, capture_output=True, text=True, timeout=300)
    except (subprocess.SubprocessError, OSError) as exc:
        raise StructureFileError(f"psfgen run failed ({driver}): {exc}") from exc
    psf = work / f"{output_prefix}.psf"
    if not psf.is_file():
        raise StructureFileError(f"psfgen did not produce {psf}")
    return psf


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate FFtk-consistent PDB (and PSF via psfgen when available)")
    parser.add_argument("mol2", type=Path)
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--topology", action="append", default=[], help="CHARMM .rtf/.str topology (repeatable); enables PSF via psfgen")
    parser.add_argument("--segid", default="MOL")
    parser.add_argument("--output-prefix", default="molecule")
    args = parser.parse_args()
    mol2 = parse_mol2(args.mol2)
    pdb_text = build_pdb(mol2)
    args.pdb.parent.mkdir(parents=True, exist_ok=True)
    args.pdb.write_text(pdb_text, encoding="utf-8")
    print(args.pdb)
    if args.topology:
        driver = psfgen_available()
        if driver is None:
            print("psfgen/VMD not found: wrote PDB only; PSF skipped")
        else:
            psf = generate_psf(args.topology, str(args.pdb), args.pdb.parent, segid=args.segid, output_prefix=args.output_prefix)
            print(psf)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
