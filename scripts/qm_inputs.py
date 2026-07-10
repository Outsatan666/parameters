from __future__ import annotations

import argparse
from pathlib import Path

from scripts.mol2 import Mol2Atom, Mol2Data, parse_mol2

# Standard element symbols, used to validate/normalise elements derived from SYBYL
# atom types. SYBYL atom types are element-based (``C.3``, ``N.ar``, ``Cl``), so the
# element is the token before the first ``.``; this set guards two-letter cases
# (Cl, Br, Na) and lets us fall back to the atom name when a type is not element-like.
ELEMENTS = frozenset(
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu "
    "Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba "
    "La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi "
    "Po At Rn Fr Ra Ac Th Pa U Np Pu".split()
)


class QMInputError(ValueError):
    """Raised when ORCA/CREST input generation receives invalid parameters."""


def _normalise(symbol: str) -> str:
    return symbol[:1].upper() + symbol[1:].lower() if symbol else symbol


def element_symbol(atom: Mol2Atom) -> str:
    """Derive a clean element symbol from a MOL2 atom.

    SYBYL atom types carry the element before the first ``.`` (``C.3`` -> ``C``,
    ``Cl`` -> ``Cl``). When the type is not element-like, fall back to the leading
    letters of the atom name (``CL1`` -> ``Cl``).
    """
    base = _normalise(atom.atom_type.split(".")[0].strip())
    if base in ELEMENTS:
        return base
    letters = "".join(ch for ch in atom.name if ch.isalpha())
    for size in (2, 1):
        candidate = _normalise(letters[:size])
        if candidate in ELEMENTS:
            return candidate
    raise QMInputError(f"Cannot derive element from atom_type={atom.atom_type!r} name={atom.name!r}")


def _coordinate_lines(mol2: Mol2Data) -> list[str]:
    return [f"{element_symbol(a):<2s} {a.x:15.6f} {a.y:15.6f} {a.z:15.6f}".rstrip() for a in mol2.atoms]


def build_orca_input(mol2: Mol2Data, *, keywords: str, charge: int = 0, multiplicity: int = 1, solvent: str | None = None) -> str:
    """Build an ORCA ``.inp`` from a parsed molecule.

    ``keywords`` is the user-supplied ``!`` line (method/basis/task); ``solvent``
    optionally appends ``CPCM(<solvent>)`` for implicit solvation. Atom order and
    element symbols come straight from the input MOL2.
    """
    if multiplicity < 1:
        raise QMInputError(f"Multiplicity must be >= 1, got {multiplicity}")
    if not mol2.atoms:
        raise QMInputError("Cannot build ORCA input for a molecule with no atoms")
    keyword_line = keywords if keywords.lstrip().startswith("!") else f"! {keywords.strip()}"
    if solvent:
        keyword_line = f"{keyword_line.rstrip()} CPCM({solvent})"
    lines = [keyword_line.rstrip(), f"* xyz {charge} {multiplicity}"]
    lines.extend(_coordinate_lines(mol2))
    lines.append("*")
    return "\n".join(lines) + "\n"


def build_crest_xyz(mol2: Mol2Data) -> str:
    """Build a standard XYZ file (atom count, comment, element x y z) for CREST/xtb."""
    if not mol2.atoms:
        raise QMInputError("Cannot build XYZ for a molecule with no atoms")
    lines = [str(len(mol2.atoms)), mol2.name]
    lines.extend(_coordinate_lines(mol2))
    return "\n".join(lines) + "\n"


def build_crest_command(xyz_name: str, *, method: str = "gfn2", solvent: str | None = None, charge: int = 0) -> str:
    """Build a reproducible CREST command line for the given XYZ file."""
    parts = ["crest", xyz_name, f"--{method}"]
    if charge:
        parts.append(f"--chrg {charge}")
    if solvent:
        parts.append(f"--alpb {solvent}")
    return " ".join(parts)


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate ORCA/CREST inputs from a MOL2 geometry")
    parser.add_argument("mol2", type=Path)
    parser.add_argument("--orca", type=Path, help="Write an ORCA .inp to this path")
    parser.add_argument("--crest-xyz", type=Path, help="Write a CREST/xtb .xyz to this path")
    parser.add_argument("--keywords", default="! R2SCAN-3c OPT freq")
    parser.add_argument("--charge", type=int, default=0)
    parser.add_argument("--multiplicity", type=int, default=1)
    parser.add_argument("--solvent")
    args = parser.parse_args()
    mol2 = parse_mol2(args.mol2)
    if args.orca:
        args.orca.parent.mkdir(parents=True, exist_ok=True)
        args.orca.write_text(build_orca_input(mol2, keywords=args.keywords, charge=args.charge, multiplicity=args.multiplicity, solvent=args.solvent), encoding="utf-8")
        print(args.orca)
    if args.crest_xyz:
        args.crest_xyz.parent.mkdir(parents=True, exist_ok=True)
        args.crest_xyz.write_text(build_crest_xyz(mol2), encoding="utf-8")
        print(args.crest_xyz)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
