from __future__ import annotations

import io
import tarfile
from pathlib import Path
from zipfile import ZipFile

import pytest
import requests

from scripts.batch_parameterize import write_review_queue
from scripts.build_archives import build_package
from scripts.inspect_package import branch_inventory, run_qc
from scripts.manifest import ManifestStore, make_job_key
from scripts.mol2 import Mol2ParseError, parse_mol2_text
from scripts.swissparam_client import SwissParamClient, SwissParamHTTPError, SwissParamProtocolError, parse_session_number, safe_extract_tar

VALID_MOL2 = """@<TRIPOS>MOLECULE
TEST
 3 2 0 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
1 C1 0.0 0.0 0.0 C.3 1 MOL 0.10
2 O1 1.0 0.0 0.0 O.3 1 MOL -0.20
3 H1 1.5 0.0 0.0 H 1 MOL 0.10
@<TRIPOS>BOND
1 1 2 1
2 2 3 1
"""


def test_parse_valid_mol2() -> None:
    data = parse_mol2_text(VALID_MOL2)
    assert data.declared_atom_count == 3
    assert data.declared_bond_count == 2
    assert data.total_charge == pytest.approx(0.0)


def test_parse_mol2_tolerates_blank_line_between_name_and_counts() -> None:
    # The review's concern: a stray blank line inside the block shifts charge_type.
    text = VALID_MOL2.replace("TEST\n 3 2 0 0 0", "TEST\n\n 3 2 0 0 0")
    data = parse_mol2_text(text)
    assert data.name == "TEST"
    assert data.declared_atom_count == 3
    assert data.charge_type == "USER_CHARGES"


def test_parse_mol2_accepts_empty_molecule_name() -> None:
    # An empty molecule-name line is legal; it must not shift the counts/charge_type fields.
    text = VALID_MOL2.replace("@<TRIPOS>MOLECULE\nTEST\n", "@<TRIPOS>MOLECULE\n\n")
    data = parse_mol2_text(text)
    assert data.name == ""
    assert data.declared_atom_count == 3
    assert data.charge_type == "USER_CHARGES"


def test_reject_malformed_mol2_without_atom_block() -> None:
    with pytest.raises(Mol2ParseError, match="ATOM"):
        parse_mol2_text(VALID_MOL2.replace("@<TRIPOS>ATOM", "@<TRIPOS>ALT_ATOM"))


def test_parse_swissparam_session_number() -> None:
    assert parse_session_number("Session number: 65720367") == "65720367"


class FakeResponse:
    def __init__(self, status_code: int, text: str = "") -> None:
        self.status_code, self.text = status_code, text
    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise requests.HTTPError(f"HTTP {self.status_code}", response=self)


class FakeSession:
    def __init__(self, statuses: list[int]) -> None:
        self.statuses, self.calls = iter(statuses), 0
    def request(self, method: str, url: str, timeout: float, **kwargs: object) -> FakeResponse:
        self.calls += 1
        return FakeResponse(next(self.statuses), "failure")


def test_detect_failed_http_response_after_bounded_retries() -> None:
    fake = FakeSession([500, 500, 500])
    client = SwissParamClient(session=fake, max_retries=2, sleeper=lambda _: None)  # type: ignore[arg-type]
    with pytest.raises(SwissParamHTTPError):
        client.check_health()
    assert fake.calls == 3


def test_submit_reuploads_full_body_on_transient_retry(tmp_path: Path) -> None:
    mol2 = tmp_path / "m.mol2"
    mol2.write_text(VALID_MOL2)
    expected = mol2.read_bytes()
    bodies: list[bytes] = []

    class RetryThenOkSession:
        def __init__(self) -> None:
            self.calls = 0

        def request(self, method: str, url: str, timeout: float, **kwargs: object) -> FakeResponse:
            self.calls += 1
            payload = kwargs["files"]["myMol2"][1]  # type: ignore[index]
            data = payload.read() if hasattr(payload, "read") else payload
            bodies.append(data)
            if self.calls == 1:
                return FakeResponse(503, "temporarily unavailable")
            return FakeResponse(200, "Session number: 424242")

    session = RetryThenOkSession()
    client = SwissParamClient(session=session, max_retries=3, sleeper=lambda _: None)  # type: ignore[arg-type]
    session_number = client.submit(mol2, approach="both")
    assert session_number == "424242"
    assert session.calls == 2
    assert bodies[0] == expected
    assert bodies[1] == expected  # retry must re-send the full MOL2, not an empty body


def test_review_queue_reports_real_state_and_branch_availability(tmp_path: Path) -> None:
    results = [{"molecule": "M", "status": "REVIEW", "state": "protonated", "match_available": "yes", "mmff_available": "no", "issues": [{"code": "MATCH_CHARGE_MISMATCH"}]}]
    path = tmp_path / "review_queue.tsv"
    write_review_queue(results, path)
    header, row = path.read_text().splitlines()[0], path.read_text().splitlines()[1]
    columns = row.split("\t")
    names = header.split("\t")
    assert columns[names.index("state")] == "protonated"
    assert columns[names.index("match_available")] == "yes"
    assert columns[names.index("mmff_available")] == "no"


def test_manifest_round_trip_and_complete_skip(tmp_path: Path) -> None:
    path = tmp_path / "manifest.json"
    key = make_job_key("abc")
    manifest = ManifestStore.load(path)
    manifest.upsert(key, molecule="A", status="TECHNICAL_PASS", session_number="123")
    assert ManifestStore.load(path).is_complete(key)


def test_generated_session_remains_resumable(tmp_path: Path) -> None:
    manifest = ManifestStore.load(tmp_path / "manifest.json")
    key = make_job_key("sha")
    manifest.upsert(key, molecule="A", status="GENERATED", session_number="999")
    assert not manifest.is_complete(key)


def test_safe_tar_extraction_rejects_traversal(tmp_path: Path) -> None:
    archive = tmp_path / "bad.tar.gz"
    with tarfile.open(archive, "w:gz") as handle:
        payload = b"owned"
        info = tarfile.TarInfo("../escape.txt")
        info.size = len(payload)
        handle.addfile(info, io.BytesIO(payload))
    with pytest.raises(SwissParamProtocolError, match="Unsafe archive"):
        safe_extract_tar(archive, tmp_path / "out")


def test_inventory_match_and_mmff_files(tmp_path: Path) -> None:
    (tmp_path / "MATCH").mkdir()
    (tmp_path / "MMFF").mkdir()
    (tmp_path / "MATCH/a.prm").write_text("x")
    (tmp_path / "MMFF/a.mol2").write_text("x")
    inventory = branch_inventory(tmp_path)
    assert inventory["MATCH"] == ["MATCH/a.prm"]
    assert inventory["MMFF"] == ["MMFF/a.mol2"]


def test_expected_charge_comparison(tmp_path: Path) -> None:
    input_path = tmp_path / "input.mol2"
    input_path.write_text(VALID_MOL2)
    (tmp_path / "MATCH").mkdir()
    (tmp_path / "MMFF").mkdir()
    (tmp_path / "MATCH/match.mol2").write_text(VALID_MOL2.replace("0.10\n2 O1", "0.20\n2 O1"))
    (tmp_path / "MMFF/mmff.mol2").write_text(VALID_MOL2)
    qc = run_qc(input_path, tmp_path, expected_charge=0)
    assert "MATCH_CHARGE_MISMATCH" in {issue["code"] for issue in qc["issues"]}
    assert qc["status"] == "REVIEW"


def test_qc_flags_match_branch_without_parameter_files(tmp_path: Path) -> None:
    input_path = tmp_path / "input.mol2"
    input_path.write_text(VALID_MOL2)
    match_dir = tmp_path / "match"
    match_dir.mkdir()
    # Looks like a MATCH branch by path, but carries no real CHARMM parameter files.
    (match_dir / "match_output.mol2").write_text(VALID_MOL2)
    qc = run_qc(input_path, tmp_path)
    codes = {issue["code"] for issue in qc["issues"]}
    assert "MATCH_PARAMETERS_MISSING" in codes
    assert qc["status"] != "TECHNICAL_PASS"


def test_qc_passes_when_match_branch_has_parameter_file(tmp_path: Path) -> None:
    input_path = tmp_path / "input.mol2"
    input_path.write_text(VALID_MOL2)
    match_dir = tmp_path / "match"
    match_dir.mkdir()
    (match_dir / "match_output.mol2").write_text(VALID_MOL2)
    (match_dir / "match.rtf").write_text("* topology\n")
    qc = run_qc(input_path, tmp_path)
    codes = {issue["code"] for issue in qc["issues"]}
    assert "MATCH_PARAMETERS_MISSING" not in codes


def test_package_zip_generation(tmp_path: Path) -> None:
    original = tmp_path / "m.mol2"
    original.write_text(VALID_MOL2)
    raw = tmp_path / "results.tar.gz"
    with tarfile.open(raw, "w:gz"):
        pass
    extracted = tmp_path / "extracted"
    (extracted / "MATCH").mkdir(parents=True)
    (extracted / "MMFF").mkdir()
    (extracted / "MATCH/m.prm").write_text("match")
    (extracted / "MMFF/m.mol2").write_text(VALID_MOL2)
    for name, content in (("comparison.tsv", "a\tb\n"), ("comparison.md", "# comparison\n"), ("qc.json", "{}"), ("issues.tsv", "code\tseverity\tmessage\tcontext\n")):
        (tmp_path / name).write_text(content)
    output = build_package("M", original_mol2=original, raw_archive=raw, extracted_root=extracted, comparison_tsv=tmp_path / "comparison.tsv", comparison_md=tmp_path / "comparison.md", qc_json=tmp_path / "qc.json", issues_tsv=tmp_path / "issues.tsv", manifest_entry={"status": "REVIEW"}, packages_dir=tmp_path / "packages", staging_root=tmp_path / "staging")
    with ZipFile(output) as archive:
        names = set(archive.namelist())
    assert {"M/input/original.mol2", "M/swissparam/raw/results.tar.gz", "M/qc/qc.json", "M/manifest.json"} <= names


def test_package_includes_derived_artifacts(tmp_path: Path) -> None:
    original = tmp_path / "m.mol2"
    original.write_text(VALID_MOL2)
    raw = tmp_path / "results.tar.gz"
    with tarfile.open(raw, "w:gz"):
        pass
    extracted = tmp_path / "extracted"
    (extracted / "MATCH").mkdir(parents=True)
    (extracted / "MATCH/m.prm").write_text("match")
    for name, content in (("comparison.tsv", "a\tb\n"), ("comparison.md", "# c\n"), ("qc.json", "{}"), ("issues.tsv", "code\n")):
        (tmp_path / name).write_text(content)
    derived = tmp_path / "derived"
    derived.mkdir()
    (derived / "m.pdb").write_text("HETATM\nEND\n")
    (derived / "m_orca.inp").write_text("! HF\n")
    output = build_package("M", original_mol2=original, raw_archive=raw, extracted_root=extracted, comparison_tsv=tmp_path / "comparison.tsv", comparison_md=tmp_path / "comparison.md", qc_json=tmp_path / "qc.json", issues_tsv=tmp_path / "issues.tsv", manifest_entry={"status": "TECHNICAL_PASS"}, packages_dir=tmp_path / "packages", staging_root=tmp_path / "staging", derived_dir=derived)
    with ZipFile(output) as archive:
        names = set(archive.namelist())
    assert {"M/derived/m.pdb", "M/derived/m_orca.inp"} <= names
