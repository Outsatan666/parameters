from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path


class Mol2ParseError(ValueError):
    """Raised when a MOL2 file violates the minimum input contract."""


@dataclass(frozen=True)
class Mol2Atom:
    atom_id: int
    name: str
    x: float
    y: float
    z: float
    atom_type: str
    subst_id: str | None
    subst_name: str | None
    charge: float | None


@dataclass(frozen=True)
class Mol2Bond:
    bond_id: int
    origin_atom_id: int
    target_atom_id: int
    bond_type: str


@dataclass(frozen=True)
class PvPorphyrinCore:
    phosphorus_atom_id: int
    nitrogen_atom_ids: tuple[int, ...]
    oxygen_atom_ids: tuple[int, ...]


@dataclass(frozen=True)
class Mol2Data:
    name: str
    declared_atom_count: int
    declared_bond_count: int
    atoms: tuple[Mol2Atom, ...]
    bonds: tuple[Mol2Bond, ...]
    charge_type: str | None

    @property
    def total_charge(self) -> float | None:
        charges = [atom.charge for atom in self.atoms]
        if not charges or any(charge is None for charge in charges):
            return None
        return float(sum(charge for charge in charges if charge is not None))

    @property
    def duplicate_atom_names(self) -> tuple[str, ...]:
        seen: set[str] = set()
        duplicates: set[str] = set()
        for atom in self.atoms:
            if atom.name in seen:
                duplicates.add(atom.name)
            seen.add(atom.name)
        return tuple(sorted(duplicates))


REQUIRED_SECTIONS = ("@<TRIPOS>MOLECULE", "@<TRIPOS>ATOM", "@<TRIPOS>BOND")


def sha256_file(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _section_bounds(lines: list[str], section: str) -> tuple[int, int]:
    starts = [index for index, line in enumerate(lines) if line.strip() == section]
    if not starts:
        raise Mol2ParseError(f"Missing required section {section}")
    start = starts[0] + 1
    end = len(lines)
    for index in range(start, len(lines)):
        if lines[index].startswith("@<TRIPOS>"):
            end = index
            break
    return start, end


def _atom_element(atom: Mol2Atom) -> str:
    return atom.atom_type.split(".", 1)[0].upper()


def detect_pv_porphyrin_core(data: Mol2Data) -> PvPorphyrinCore | None:
    atoms_by_id = {atom.atom_id: atom for atom in data.atoms}
    neighbors: dict[int, set[int]] = {atom.atom_id: set() for atom in data.atoms}
    for bond in data.bonds:
        neighbors[bond.origin_atom_id].add(bond.target_atom_id)
        neighbors[bond.target_atom_id].add(bond.origin_atom_id)

    for atom in data.atoms:
        if _atom_element(atom) != "P":
            continue
        neighbor_ids = neighbors[atom.atom_id]
        if len(neighbor_ids) != 6:
            continue
        nitrogen_ids = tuple(sorted(atom_id for atom_id in neighbor_ids if _atom_element(atoms_by_id[atom_id]) == "N"))
        oxygen_ids = tuple(sorted(atom_id for atom_id in neighbor_ids if _atom_element(atoms_by_id[atom_id]) == "O"))
        if len(nitrogen_ids) == 4 and len(oxygen_ids) == 2:
            return PvPorphyrinCore(atom.atom_id, nitrogen_ids, oxygen_ids)
    return None


def parse_mol2_text(text: str) -> Mol2Data:
    lines = text.splitlines()
    present = {line.strip() for line in lines if line.startswith("@<TRIPOS>")}
    missing = [section for section in REQUIRED_SECTIONS if section not in present]
    if missing:
        raise Mol2ParseError(f"Missing required MOL2 sections: {', '.join(missing)}")
    molecule_start, molecule_end = _section_bounds(lines, "@<TRIPOS>MOLECULE")
    molecule_block = [line.strip() for line in lines[molecule_start:molecule_end]]
    if len(molecule_block) < 2:
        raise Mol2ParseError("MOLECULE block does not contain name and counts")
    name = molecule_block[0]
    count_fields = molecule_block[1].split()
    if len(count_fields) < 2:
        raise Mol2ParseError("MOLECULE counts line is malformed")
    try:
        declared_atom_count = int(count_fields[0])
        declared_bond_count = int(count_fields[1])
    except ValueError as exc:
        raise Mol2ParseError("MOLECULE atom/bond counts are not integers") from exc
    charge_type = molecule_block[3] if len(molecule_block) >= 4 and molecule_block[3] else None
    atom_start, atom_end = _section_bounds(lines, "@<TRIPOS>ATOM")
    atoms: list[Mol2Atom] = []
    for line_number, line in enumerate(lines[atom_start:atom_end], start=atom_start + 1):
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) < 6:
            raise Mol2ParseError(f"Malformed ATOM line {line_number}: {line!r}")
        try:
            atom_id = int(fields[0])
            x, y, z = map(float, fields[2:5])
            charge = float(fields[8]) if len(fields) >= 9 else None
        except ValueError as exc:
            raise Mol2ParseError(f"Non-numeric ATOM field on line {line_number}") from exc
        atoms.append(Mol2Atom(atom_id, fields[1], x, y, z, fields[5], fields[6] if len(fields) >= 7 else None, fields[7] if len(fields) >= 8 else None, charge))
    bond_start, bond_end = _section_bounds(lines, "@<TRIPOS>BOND")
    bonds: list[Mol2Bond] = []
    for line_number, line in enumerate(lines[bond_start:bond_end], start=bond_start + 1):
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) < 4:
            raise Mol2ParseError(f"Malformed BOND line {line_number}: {line!r}")
        try:
            bond_id, origin, target = map(int, fields[:3])
        except ValueError as exc:
            raise Mol2ParseError(f"Non-numeric BOND field on line {line_number}") from exc
        bonds.append(Mol2Bond(bond_id, origin, target, fields[3]))
    if len(atoms) != declared_atom_count:
        raise Mol2ParseError(f"Declared atom count {declared_atom_count} != parsed atom count {len(atoms)}")
    if len(bonds) != declared_bond_count:
        raise Mol2ParseError(f"Declared bond count {declared_bond_count} != parsed bond count {len(bonds)}")
    atom_ids = [atom.atom_id for atom in atoms]
    if len(atom_ids) != len(set(atom_ids)):
        raise Mol2ParseError("Duplicate atom IDs detected")
    return Mol2Data(name, declared_atom_count, declared_bond_count, tuple(atoms), tuple(bonds), charge_type)


def parse_mol2(path: str | Path) -> Mol2Data:
    return parse_mol2_text(Path(path).read_text(encoding="utf-8", errors="replace"))
