from __future__ import annotations

import os
from pathlib import Path

import pytest

from scripts.swissparam_client import SwissParamClient


@pytest.mark.integration
def test_swissparam_health_endpoint() -> None:
    if os.environ.get("RUN_SWISSPARAM_INTEGRATION") != "1":
        pytest.skip("Set RUN_SWISSPARAM_INTEGRATION=1 for explicit external integration test")
    assert "Hello World!" in SwissParamClient().check_health()


@pytest.mark.integration
def test_swissparam_single_molecule_submission(tmp_path: Path) -> None:
    if os.environ.get("RUN_SWISSPARAM_INTEGRATION") != "1":
        pytest.skip("Set RUN_SWISSPARAM_INTEGRATION=1 for explicit external integration test")
    mol2_env = os.environ.get("SWISSPARAM_TEST_MOL2")
    if not mol2_env:
        pytest.skip("Set SWISSPARAM_TEST_MOL2 to an explicit test MOL2 path")
    client = SwissParamClient()
    session = client.submit(Path(mol2_env), approach="both")
    result = client.poll(session, poll_interval_s=15, max_total_wait_s=7200)
    assert result.state == "finished"
    assert client.retrieve(session, tmp_path / "results.tar.gz").is_file()
