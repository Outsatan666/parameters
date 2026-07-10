from __future__ import annotations

import subprocess

import pytest

import scripts.structure_files as sf
from scripts.mol2 import parse_mol2_text
from scripts.structure_files import (
    StructureFileError,
    build_pdb,
    build_psfgen_script,
    generate_psf,
    numbering_matches,
)

VALID_MOL2 = """@<TRIPOS>MOLECULE
LIG
 3 2 0 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
1 C1 0.0000 0.0000 0.0000 C.3 1 LIG 0.10
2 O1 1.2000 0.0000 0.0000 O.3 1 LIG -0.20
3 H1 1.5000 0.9000 0.0000 H 1 LIG 0.10
@<TRIPOS>BOND
1 1 2 1
2 2 3 1
"""


def test_pdb_has_one_record_per_atom_ending_with_end() -> None:
    mol2 = parse_mol2_text(VALID_MOL2)
    lines = build_pdb(mol2).splitlines()
    atom_lines = [ln for ln in lines if ln.startswith(("ATOM", "HETATM"))]
    assert len(atom_lines) == 3
    assert lines[-1] == "END"


def test_pdb_serial_and_name_match_mol2_atom_numbering() -> None:
    mol2 = parse_mol2_text(VALID_MOL2)
    atom_lines = [ln for ln in build_pdb(mol2).splitlines() if ln.startswith(("ATOM", "HETATM"))]
    for atom, line in zip(mol2.atoms, atom_lines, strict=True):
        serial = int(line[6:11])
        name = line[12:16].strip()
        assert serial == atom.atom_id
        assert name == atom.name


def test_pdb_columns_carry_element_and_coordinates() -> None:
    mol2 = parse_mol2_text(VALID_MOL2)
    first = [ln for ln in build_pdb(mol2).splitlines() if ln.startswith(("ATOM", "HETATM"))][0]
    assert first[76:78].strip() == "C"
    assert float(first[30:38]) == pytest.approx(0.0)
    assert float(first[38:46]) == pytest.approx(0.0)


def test_pdb_uses_mol2_residue_name_by_default() -> None:
    mol2 = parse_mol2_text(VALID_MOL2)
    first = [ln for ln in build_pdb(mol2).splitlines() if ln.startswith(("ATOM", "HETATM"))][0]
    assert first[17:20].strip() == "LIG"


def test_build_pdb_rejects_atom_name_longer_than_four_chars() -> None:
    # PDB has only 4 name columns; silently truncating would break the mol2<->pdb name identity.
    mol2 = parse_mol2_text(VALID_MOL2.replace(" C1 ", " CARBON "))
    with pytest.raises(StructureFileError):
        build_pdb(mol2)


def test_numbering_matches_is_true_for_generated_pdb() -> None:
    mol2 = parse_mol2_text(VALID_MOL2)
    assert numbering_matches(mol2, build_pdb(mol2)) is True


def test_numbering_matches_detects_reordered_atoms() -> None:
    mol2 = parse_mol2_text(VALID_MOL2)
    pdb = build_pdb(mol2)
    lines = pdb.splitlines()
    atom_idx = [i for i, ln in enumerate(lines) if ln.startswith(("ATOM", "HETATM"))]
    lines[atom_idx[0]], lines[atom_idx[1]] = lines[atom_idx[1]], lines[atom_idx[0]]
    assert numbering_matches(mol2, "\n".join(lines)) is False


def test_psfgen_script_references_topology_pdb_and_segid() -> None:
    script = build_psfgen_script(["match.rtf"], "molecule.pdb", segid="LIG", output_prefix="molecule")
    assert "topology match.rtf" in script
    assert "coordpdb molecule.pdb LIG" in script
    assert "writepsf molecule.psf" in script
    assert "writepdb molecule.pdb" in script


def test_psfgen_script_requires_at_least_one_topology() -> None:
    with pytest.raises(StructureFileError):
        build_psfgen_script([], "molecule.pdb", segid="LIG")


def test_generate_psf_wraps_subprocess_failure(monkeypatch: pytest.MonkeyPatch, tmp_path) -> None:
    monkeypatch.setattr(sf, "psfgen_available", lambda: "vmd")

    def boom(*args: object, **kwargs: object) -> object:
        raise subprocess.CalledProcessError(1, "vmd", stderr="psfgen failed")

    monkeypatch.setattr(sf.subprocess, "run", boom)
    topology = tmp_path / "match.rtf"
    topology.write_text("* topology\n")
    pdb = tmp_path / "m.pdb"
    pdb.write_text("HETATM\nEND\n")
    with pytest.raises(StructureFileError):
        generate_psf([str(topology)], str(pdb), tmp_path, output_prefix="m")


def test_generate_psf_wraps_timeout(monkeypatch: pytest.MonkeyPatch, tmp_path) -> None:
    monkeypatch.setattr(sf, "psfgen_available", lambda: "psfgen")

    def slow(*args: object, **kwargs: object) -> object:
        raise subprocess.TimeoutExpired("psfgen", 300)

    monkeypatch.setattr(sf.subprocess, "run", slow)
    topology = tmp_path / "match.str"
    topology.write_text("* topology\n")
    pdb = tmp_path / "m.pdb"
    pdb.write_text("HETATM\nEND\n")
    with pytest.raises(StructureFileError):
        generate_psf([str(topology)], str(pdb), tmp_path, output_prefix="m")
