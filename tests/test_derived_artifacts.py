from __future__ import annotations

from pathlib import Path

from scripts.derived_artifacts import write_derived_artifacts

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


def _mol2(tmp_path: Path) -> Path:
    p = tmp_path / "LIG.mol2"
    p.write_text(VALID_MOL2)
    return p


def test_writes_pdb_orca_and_crest_artifacts(tmp_path: Path) -> None:
    mol2 = _mol2(tmp_path)
    out = tmp_path / "derived"
    written = write_derived_artifacts(mol2, out, orca_keywords="! R2SCAN-3c OPT freq")
    assert (out / "LIG.pdb").is_file()
    assert (out / "LIG_orca.inp").is_file()
    assert (out / "LIG.xyz").is_file()
    assert "pdb" in written and "orca" in written and "crest_xyz" in written


def test_orca_artifact_uses_requested_keywords_and_charge(tmp_path: Path) -> None:
    mol2 = _mol2(tmp_path)
    out = tmp_path / "derived"
    write_derived_artifacts(mol2, out, orca_keywords="! B3LYP def2-SVP", charge=-1, multiplicity=1, solvent="Water")
    inp = (out / "LIG_orca.inp").read_text().splitlines()
    assert inp[0] == "! B3LYP def2-SVP CPCM(Water)"
    assert inp[1] == "* xyz -1 1"


def test_pdb_numbering_matches_input_mol2(tmp_path: Path) -> None:
    from scripts.mol2 import parse_mol2
    from scripts.structure_files import numbering_matches

    mol2 = _mol2(tmp_path)
    out = tmp_path / "derived"
    write_derived_artifacts(mol2, out)
    assert numbering_matches(parse_mol2(mol2), (out / "LIG.pdb").read_text()) is True


def test_writes_runnable_crest_command_matching_solvent_and_charge(tmp_path: Path) -> None:
    mol2 = _mol2(tmp_path)
    out = tmp_path / "derived"
    written = write_derived_artifacts(mol2, out, solvent="Water", charge=-1)
    command = (out / "LIG_crest.sh").read_text()
    assert "crest LIG.xyz" in command
    assert "--alpb water" in command
    assert "--chrg -1" in command
    assert "crest_cmd" in written


def test_psf_skipped_gracefully_when_psfgen_absent(tmp_path: Path) -> None:
    # No psfgen/VMD on PATH in CI/dev: PSF must be skipped, not crash.
    mol2 = _mol2(tmp_path)
    out = tmp_path / "derived"
    written = write_derived_artifacts(mol2, out, topology_paths=["nonexistent.rtf"])
    assert "psf" not in written or written["psf"] is None
