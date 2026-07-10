from __future__ import annotations

import json
import os
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any

PIPELINE_VERSION = "0.1.0"
SUCCESS_STATUSES = {"TECHNICAL_PASS", "SECOND_LEVEL_COMPLETE"}


def make_job_key(input_sha256: str, generator: str = "swissparam", approach: str = "both", pipeline_version: str = PIPELINE_VERSION) -> str:
    payload = "\n".join((input_sha256, generator, approach, pipeline_version)).encode("utf-8")
    return sha256(payload).hexdigest()


@dataclass
class ManifestStore:
    path: Path
    data: dict[str, Any]

    @classmethod
    def load(cls, path: str | Path) -> "ManifestStore":
        manifest_path = Path(path)
        if manifest_path.exists():
            data = json.loads(manifest_path.read_text(encoding="utf-8"))
        else:
            data = {"schema_version": 1, "pipeline_version": PIPELINE_VERSION, "entries": {}}
        data.setdefault("entries", {})
        return cls(manifest_path, data)

    def get(self, job_key: str) -> dict[str, Any] | None:
        entry = self.data["entries"].get(job_key)
        return dict(entry) if entry is not None else None

    def is_complete(self, job_key: str) -> bool:
        entry = self.data["entries"].get(job_key)
        return bool(entry and entry.get("status") in SUCCESS_STATUSES)

    def upsert(self, job_key: str, **fields: Any) -> dict[str, Any]:
        entry = self.data["entries"].setdefault(job_key, {})
        entry.update(fields)
        self.save()
        return dict(entry)

    def save(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = self.path.with_suffix(self.path.suffix + ".tmp")
        tmp_path.write_text(json.dumps(self.data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(tmp_path, self.path)
