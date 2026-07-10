from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path
from typing import Protocol

OUTPUTS_ROOT = "outputs"


class UploadError(RuntimeError):
    """Raised when publishing results to a destination fails or is misconfigured."""


class Uploader(Protocol):
    def publish(self, molecule: str, files: list[Path]) -> list[str]:
        """Publish ``files`` for ``molecule`` and return destination identifiers."""


class LocalDirectoryUploader:
    """Publish into ``<root>/outputs/<molecule>/``.

    Works today and is the credential-free way to reach Google Drive: point ``root``
    at a Drive-mounted folder or an rclone/Insync staging directory. This code never
    handles Drive credentials.
    """

    def __init__(self, root: str | Path) -> None:
        self.root = Path(root)

    def publish(self, molecule: str, files: list[Path]) -> list[Path]:
        target_dir = self.root / OUTPUTS_ROOT / molecule
        target_dir.mkdir(parents=True, exist_ok=True)
        published: list[Path] = []
        for file in files:
            source = Path(file)
            if not source.is_file():
                raise UploadError(f"Cannot publish missing file: {source}")
            destination = target_dir / source.name
            shutil.copy2(source, destination)
            published.append(destination)
        return published


class RcloneUploader:
    """Publish via ``rclone`` to ``<remote>/outputs/<molecule>/``.

    Credentials live in the user's own rclone configuration, never here. Requires
    ``rclone`` on PATH and a remote the user has already configured (e.g. ``gdrive:``).
    """

    def __init__(self, remote: str) -> None:
        self.remote = remote.rstrip("/")

    def publish(self, molecule: str, files: list[Path]) -> list[str]:
        rclone = shutil.which("rclone")
        if rclone is None:
            raise UploadError("rclone not found on PATH; configure rclone or use a local: destination")
        published: list[str] = []
        for file in files:
            source = Path(file)
            if not source.is_file():
                raise UploadError(f"Cannot publish missing file: {source}")
            destination = f"{self.remote}/{OUTPUTS_ROOT}/{molecule}/{source.name}"
            try:
                subprocess.run([rclone, "copyto", str(source), destination], check=True, capture_output=True, text=True, timeout=600)
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
                raise UploadError(f"rclone failed to publish {source} to {destination}: {exc}") from exc
            published.append(destination)
        return published


def resolve_uploader(spec: str | None) -> Uploader | None:
    """Build an uploader from a spec, or None when publishing is not configured.

    Specs: ``local:<path>`` (copy into a directory) or ``rclone:<remote>`` (e.g.
    ``rclone:gdrive:CHARMM``). An empty/None spec disables publishing.
    """
    if not spec:
        return None
    scheme, _, target = spec.partition(":")
    if scheme == "local":
        return LocalDirectoryUploader(target)
    if scheme == "rclone":
        if not target:
            raise UploadError("rclone spec requires a remote, e.g. rclone:gdrive:CHARMM")
        return RcloneUploader(target)
    raise UploadError(f"Unknown upload scheme {scheme!r}; use local:<path> or rclone:<remote>")


def main() -> int:
    parser = argparse.ArgumentParser(description="Publish result files for a molecule to a destination")
    parser.add_argument("molecule")
    parser.add_argument("files", nargs="+", type=Path)
    parser.add_argument("--destination", required=True, help="local:<path> or rclone:<remote>")
    args = parser.parse_args()
    uploader = resolve_uploader(args.destination)
    if uploader is None:
        print("No destination configured; nothing published")
        return 0
    for published in uploader.publish(args.molecule, args.files):
        print(published)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
