from __future__ import annotations

from pathlib import Path

import pytest

import scripts.upload as up
from scripts.upload import (
    LocalDirectoryUploader,
    RcloneUploader,
    UploadError,
    resolve_uploader,
)


def test_local_uploader_copies_files_under_outputs_molecule(tmp_path: Path) -> None:
    package = tmp_path / "M_CHARMM_STARTING_PARAMETERS.zip"
    package.write_bytes(b"zip-bytes")
    dest = tmp_path / "drive"
    uploader = LocalDirectoryUploader(dest)
    published = uploader.publish("M", [package])
    target = dest / "outputs" / "M" / package.name
    assert target.is_file()
    assert target.read_bytes() == b"zip-bytes"
    assert published == [target]


def test_local_uploader_creates_nested_destination(tmp_path: Path) -> None:
    package = tmp_path / "p.zip"
    package.write_bytes(b"x")
    dest = tmp_path / "a" / "b" / "drive"
    LocalDirectoryUploader(dest).publish("LIG", [package])
    assert (dest / "outputs" / "LIG" / "p.zip").is_file()


def test_resolve_uploader_returns_none_when_unset() -> None:
    assert resolve_uploader(None) is None
    assert resolve_uploader("") is None


def test_resolve_uploader_local_spec(tmp_path: Path) -> None:
    uploader = resolve_uploader(f"local:{tmp_path / 'out'}")
    assert isinstance(uploader, LocalDirectoryUploader)


def test_resolve_uploader_rejects_unknown_scheme() -> None:
    with pytest.raises(UploadError):
        resolve_uploader("ftp:whatever")


def test_resolve_uploader_rclone_spec_builds_rclone_uploader() -> None:
    assert isinstance(resolve_uploader("rclone:gdrive:CHARMM"), RcloneUploader)


def test_resolve_uploader_rclone_requires_remote() -> None:
    with pytest.raises(UploadError):
        resolve_uploader("rclone:")


def test_rclone_uploader_builds_copyto_command_and_destination(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(up.shutil, "which", lambda name: "/usr/bin/rclone" if name == "rclone" else None)
    monkeypatch.setattr(up.subprocess, "run", lambda cmd, **kwargs: calls.append(cmd))
    package = tmp_path / "p.zip"
    package.write_bytes(b"x")
    published = RcloneUploader("gdrive:CHARMM").publish("LIG", [package])
    assert calls[0][:3] == ["/usr/bin/rclone", "copyto", str(package)]
    assert calls[0][3] == "gdrive:CHARMM/outputs/LIG/p.zip"
    assert published == ["gdrive:CHARMM/outputs/LIG/p.zip"]


def test_rclone_uploader_errors_when_rclone_absent(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setattr(up.shutil, "which", lambda name: None)
    package = tmp_path / "p.zip"
    package.write_bytes(b"x")
    with pytest.raises(UploadError):
        RcloneUploader("gdrive:X").publish("LIG", [package])


def test_rclone_uploader_wraps_subprocess_failure(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    import subprocess

    monkeypatch.setattr(up.shutil, "which", lambda name: "/usr/bin/rclone")

    def boom(cmd: list[str], **kwargs: object) -> object:
        raise subprocess.CalledProcessError(1, cmd)

    monkeypatch.setattr(up.subprocess, "run", boom)
    package = tmp_path / "p.zip"
    package.write_bytes(b"x")
    with pytest.raises(UploadError):
        RcloneUploader("gdrive:X").publish("LIG", [package])
