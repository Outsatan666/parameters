from __future__ import annotations

import pytest

from scripts.mol2 import parse_mol2_text
from scripts.qm_inputs import (
    QMInputError,
    build_crest_xyz,
    build_orca_input,
    element_symbol,
)

VALID_MOL2 = """@<TRIPOS>MOLECULE
TEST
 3 2 0 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
1 C1 0.0000 0.0000 0.0000 C.3 1 MOL 0.10
2 O1 1.2000 0.0000 0.0000 O.3 1 MOL -0.20
3 H1 1.5000 0.9000 0.0000 H 1 MOL 0.10
@<TRIPOS>BOND
1 1 2 1
2 2 3 1
"""

HALOGEN_MOL2 = """@<TRIPOS>MOLECULE
HAL
 2 1 0 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
1 C1 0.0 0.0 0.0 C.3 1 MOL 0.10
2 CL1 1.7 0.0 0.0 Cl 1 MOL -0.10
@<TRIPOS>BOND
1 1 2 1
"""


def _atoms(text: str):
    return parse_mol2_text(text).atoms


def test_element_symbol_from_sybyl_types() -> None:
    atoms = _atoms(VALID_MOL2)
    assert [element_symbol(a) for a in atoms] == ["C", "O", "H"]


def test_element_symbol_handles_two_letter_halogen() -> None:
    atoms = _atoms(HALOGEN_MOL2)
    assert element_symbol(atoms[1]) == "Cl"


def test_element_symbol_raises_for_non_element_atom() -> None:
    # Dummy atoms / lone pairs are not real elements: fail loudly rather than emit a bogus atom.
    dummy = """@<TRIPOS>MOLECULE
D
 1 0 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
1 XX 0.0 0.0 0.0 Du 1 MOL 0.0
@<TRIPOS>BOND
"""
    with pytest.raises(QMInputError):
        element_symbol(_atoms(dummy)[0])


def test_orca_input_default_charge_multiplicity_and_coords() -> None:
    mol2 = parse_mol2_text(VALID_MOL2)
    text = build_orca_input(mol2, keywords="! R2SCAN-3c OPT freq")
    lines = text.splitlines()
    assert lines[0] == "! R2SCAN-3c OPT freq"
    assert lines[1] == "* xyz 0 1"
    assert lines[2].split() == ["C", "0.000000", "0.000000", "0.000000"]
    assert lines[3].split()[0] == "O"
    assert lines[4].split()[0] == "H"
    assert lines[5] == "*"
    # one coordinate line per atom, plus keyword/charge/closing lines
    assert sum(1 for line in lines if line and line[0].isalpha() and not line.startswith("!")) == 3


def test_orca_input_overrides_charge_and_multiplicity() -> None:
    mol2 = parse_mol2_text(VALID_MOL2)
    text = build_orca_input(mol2, keywords="! B3LYP def2-SVP", charge=-1, multiplicity=2)
    assert "* xyz -1 2" in text.splitlines()


def test_orca_input_appends_water_solvation_when_requested() -> None:
    mol2 = parse_mol2_text(VALID_MOL2)
    text = build_orca_input(mol2, keywords="! R2SCAN-3c OPT", solvent="Water")
    assert text.splitlines()[0] == "! R2SCAN-3c OPT CPCM(Water)"


def test_orca_input_rejects_nonpositive_multiplicity() -> None:
    mol2 = parse_mol2_text(VALID_MOL2)
    with pytest.raises(QMInputError):
        build_orca_input(mol2, keywords="! HF", multiplicity=0)


def test_crest_xyz_is_standard_xyz_with_count_and_comment() -> None:
    mol2 = parse_mol2_text(VALID_MOL2)
    text = build_crest_xyz(mol2)
    lines = text.splitlines()
    assert lines[0] == "3"
    assert lines[1] == "TEST"
    assert lines[2].split() == ["C", "0.000000", "0.000000", "0.000000"]
    assert len(lines) == 5  # count + comment + 3 atoms


def test_orca_and_crest_preserve_input_atom_order() -> None:
    mol2 = parse_mol2_text(VALID_MOL2)
    orca_elements = [ln.split()[0] for ln in build_orca_input(mol2, keywords="! HF").splitlines()[2:5]]
    crest_elements = [ln.split()[0] for ln in build_crest_xyz(mol2).splitlines()[2:5]]
    assert orca_elements == ["C", "O", "H"]
    assert crest_elements == ["C", "O", "H"]
