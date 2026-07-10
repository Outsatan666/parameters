from __future__ import annotations

import io
import tarfile
from pathlib import Path
from zipfile import ZipFile

import pytest
import requests

from scripts.batch_parameterize import process_one
from scripts.build_archives import build_package
from scripts.inspect_package import branch_inventory, run_qc
from scripts.manifest import ManifestStore, make_job_key
from scripts.mol2 import Mol2ParseError, detect_pv_porphyrin_core, parse_mol2_text
from scripts.swissparam_client import SwissParamClient, SwissParamHTTPError, SwissParamProtocolError, classify_status, parse_session_number, safe_extract_tar, summarize_failure

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

PV_PORPHYRIN_MOL2 = """@<TRIPOS>MOLECULE
PV_CORE
 7 6 0 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
1 P 0.0 0.0 0.0 P.3 1 MOL 0.20
2 N1 1.0 0.0 0.0 N.2 1 MOL -0.05
3 N2 -1.0 0.0 0.0 N.2 1 MOL -0.05
4 N3 0.0 1.0 0.0 N.ar 1 MOL -0.05
5 N4 0.0 -1.0 0.0 N.ar 1 MOL -0.05
6 O1 0.0 0.0 1.0 O.3 1 MOL 0.00
7 O2 0.0 0.0 -1.0 O.3 1 MOL 0.00
@<TRIPOS>BOND
1 1 2 1
2 1 3 1
3 1 4 1
4 1 5 1
5 1 6 1
6 1 7 1
"""

PYROH2_SWISSPARAM_FAILURE = """
Your parameters:
    Query: mol2
    Input: PyrOH2.mol2
    Approach: both
Could not retrieve alternative tautomers and protonation states.
ERROR: Can't use an undefined value as an ARRAY reference
ERROR. Problem in CHARMM run. Check tmp.out.
ERROR PyrOH2.mol2 COULD NOT BE DONE (8)
ERROR_pka, Error occurred while processing mol2 file: Python argument types in
"""


def test_parse_valid_mol2() -> None:
    data = parse_mol2_text(VALID_MOL2)
    assert data.declared_atom_count == 3
    assert data.declared_bond_count == 2
    assert data.total_charge == pytest.approx(0.0)


def test_reject_malformed_mol2_without_atom_block() -> None:
    with pytest.raises(Mol2ParseError, match="ATOM"):
        parse_mol2_text(VALID_MOL2.replace("@<TRIPOS>ATOM", "@<TRIPOS>ALT_ATOM"))


def test_detect_pv_porphyrin_core_by_p_n4_o2_coordination() -> None:
    core = detect_pv_porphyrin_core(parse_mol2_text(PV_PORPHYRIN_MOL2))
    assert core is not None
    assert core.phosphorus_atom_id == 1
    assert core.nitrogen_atom_ids == (2, 3, 4, 5)
    assert core.oxygen_atom_ids == (6, 7)


def test_ordinary_mol2_does_not_trigger_pv_porphyrin_route() -> None:
    assert detect_pv_porphyrin_core(parse_mol2_text(VALID_MOL2)) is None


def test_parse_swissparam_session_number() -> None:
    assert parse_session_number("Session number: 65720367") == "65720367"


def test_classify_actual_pyroh2_swissparam_failure() -> None:
    assert classify_status(PYROH2_SWISSPARAM_FAILURE) == "failed"


def test_terminal_failure_takes_precedence_over_running_text() -> None:
    text = "Calculation currently running.\nERROR PyrOH2.mol2 COULD NOT BE DONE (8)"
    assert classify_status(text) == "failed"


def test_summarize_actual_pyroh2_failure_preserves_backend_reason() -> None:
    summary = summarize_failure(PYROH2_SWISSPARAM_FAILURE)
    assert "Problem in CHARMM run" in summary
    assert "COULD NOT BE DONE (8)" in summary
    assert "exceeded max wait" not in summary


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


class StaticTextSession:
    def __init__(self, text: str) -> None:
        self.text, self.calls = text, 0
    def request(self, method: str, url: str, timeout: float, **kwargs: object) -> FakeResponse:
        self.calls += 1
        return FakeResponse(200, self.text)


class RejectingSwissParamClient:
    def check_health(self) -> str:
        raise AssertionError("specialized P(V)-porphyrin route must not call SwissParam health")
    def submit(self, *args: object, **kwargs: object) -> str:
        raise AssertionError("specialized P(V)-porphyrin route must not submit to SwissParam")
    def poll(self, *args: object, **kwargs: object) -> object:
        raise AssertionError("specialized P(V)-porphyrin route must not poll SwissParam")


def test_detect_failed_http_response_after_bounded_retries() -> None:
    fake = FakeSession([500, 500, 500])
    client = SwissParamClient(session=fake, max_retries=2, sleeper=lambda _: None)  # type: ignore[arg-type]
    with pytest.raises(SwissParamHTTPError):
        client.check_health()
    assert fake.calls == 3


def test_poll_stops_on_actual_pyroh2_failure_without_sleeping() -> None:
    fake = StaticTextSession(PYROH2_SWISSPARAM_FAILURE)
    sleeps: list[float] = []
    client = SwissParamClient(session=fake, sleeper=sleeps.append)  # type: ignore[arg-type]
    result = client.poll("74281421", poll_interval_s=15, max_total_wait_s=7200)
    assert result.state == "failed"
    assert fake.calls == 1
    assert sleeps == []


def test_pv_porphyrin_routes_to_custom_ff_before_swissparam(tmp_path: Path) -> None:
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    mol2 = input_dir / "PyrOH2.mol2"
    mol2.write_text(PV_PORPHYRIN_MOL2)
    manifest = ManifestStore.load(tmp_path / "manifests/manifest.json")

    result = process_one(
        mol2,
        client=RejectingSwissParamClient(),  # type: ignore[arg-type]
        manifest=manifest,
        state_row={"expected_charge": "1"},
        approach="both",
        poll_interval_s=15,
        max_total_wait_s=7200,
        force=False,
        workspace=tmp_path,
    )

    assert result["status"] == "REVIEW"
    assert result["state"] == "SPECIALIZED_PARAMETERIZATION_REQUIRED"
    assert result["action"] == "ROUTE_SPECIALIZED_PARAMETERIZATION"
    assert result["review_reason"] == "PV_PORPHYRIN_CUSTOM_FF_REQUIRED"
    assert result["specialized_route"]["phosphorus_atom_id"] == 1
    assert result["specialized_route"]["neighbor_element_counts"] == {"N": 4, "O": 2}
    issue_codes = {issue["code"] for issue in result["issues"]}
    assert "PV_PORPHYRIN_CUSTOM_FF_REQUIRED" in issue_codes
    assert "INPUT_CHARGE_MISMATCH" in issue_codes


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
