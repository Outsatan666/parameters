from __future__ import annotations

import io
import tarfile
from pathlib import Path
from zipfile import ZipFile

import pytest
import requests

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
