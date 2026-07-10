from __future__ import annotations

import argparse
import re
import tarfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import requests

DEFAULT_BASE_URL = "https://www.swissparam.ch:8443"
SESSION_RE = re.compile(r"Session\s+number\s*:\s*(\d+)", re.IGNORECASE)
FAILURE_MARKERS = (
    "calculation failed",
    "has failed",
    "fatal error",
    "could not be done",
    "error. problem in charmm run",
    "error_pka",
    "error occurred while processing mol2 file",
)


class SwissParamError(RuntimeError):
    pass


class SwissParamHTTPError(SwissParamError):
    pass


class SwissParamProtocolError(SwissParamError):
    pass


@dataclass(frozen=True)
class StatusResult:
    state: str
    raw_text: str


def parse_session_number(text: str) -> str:
    match = SESSION_RE.search(text)
    if not match:
        raise SwissParamProtocolError("SwissParam response did not contain a Session number")
    return match.group(1)


def classify_status(text: str) -> str:
    lowered = text.lower()
    if "calculation is finished" in lowered:
        return "finished"
    if any(token in lowered for token in FAILURE_MARKERS):
        return "failed"
    if "currently running" in lowered:
        return "running"
    if "in the queue" in lowered or "pending" in lowered:
        return "queued"
    return "unknown"


def safe_extract_tar(archive_path: str | Path, destination: str | Path) -> list[Path]:
    archive = Path(archive_path)
    destination_path = Path(destination).resolve()
    destination_path.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:*") as handle:
        for member in handle.getmembers():
            target = (destination_path / member.name).resolve()
            try:
                target.relative_to(destination_path)
            except ValueError as exc:
                raise SwissParamProtocolError(f"Unsafe archive member path detected: {member.name}") from exc
        handle.extractall(destination_path, filter="data")
    return [path for path in destination_path.rglob("*") if path.is_file()]


class SwissParamClient:
    def __init__(self, base_url: str = DEFAULT_BASE_URL, *, timeout_s: float = 60.0, max_retries: int = 4, retry_backoff_s: float = 2.0, session: requests.Session | None = None, sleeper: Callable[[float], None] = time.sleep) -> None:
        self.base_url = base_url.rstrip("/")
        self.timeout_s = timeout_s
        self.max_retries = max_retries
        self.retry_backoff_s = retry_backoff_s
        self.session = session or requests.Session()
        self.sleeper = sleeper

    def _request(self, method: str, url: str, **kwargs: object) -> requests.Response:
        last_error: Exception | None = None
        for attempt in range(self.max_retries + 1):
            try:
                response = self.session.request(method, url, timeout=self.timeout_s, **kwargs)
                if response.status_code >= 500 or response.status_code in {408, 429}:
                    raise requests.HTTPError(f"Transient HTTP {response.status_code}", response=response)
                response.raise_for_status()
                return response
            except (requests.RequestException, OSError) as exc:
                last_error = exc
                if attempt >= self.max_retries:
                    break
                self.sleeper(self.retry_backoff_s * (2**attempt))
        raise SwissParamHTTPError(f"Request failed after {self.max_retries + 1} attempts: {last_error}") from last_error

    def check_health(self) -> str:
        text = self._request("GET", f"{self.base_url}/").text.strip()
        if "Hello World!" not in text:
            raise SwissParamProtocolError(f"Unexpected SwissParam health response: {text[:200]!r}")
        return text

    def submit(self, mol2_path: str | Path, *, approach: str = "both", raw_response_path: str | Path | None = None) -> str:
        if approach not in {"both", "mmff-based", "match"}:
            raise ValueError(f"Unsupported SwissParam approach: {approach}")
        path = Path(mol2_path)
        with path.open("rb") as handle:
            response = self._request("POST", f"{self.base_url}/startparam", params={"approach": approach}, files={"myMol2": (path.name, handle, "chemical/x-mol2")})
        if raw_response_path is not None:
            raw = Path(raw_response_path)
            raw.parent.mkdir(parents=True, exist_ok=True)
            raw.write_text(response.text, encoding="utf-8")
        return parse_session_number(response.text)

    def check_status(self, session_number: str, *, raw_response_path: str | Path | None = None) -> StatusResult:
        response = self._request("GET", f"{self.base_url}/checksession", params={"sessionNumber": session_number})
        if raw_response_path is not None:
            raw = Path(raw_response_path)
            raw.parent.mkdir(parents=True, exist_ok=True)
            raw.write_text(response.text, encoding="utf-8")
        return StatusResult(classify_status(response.text), response.text)

    def poll(self, session_number: str, *, poll_interval_s: float = 15.0, max_total_wait_s: float = 7200.0, raw_status_dir: str | Path | None = None) -> StatusResult:
        deadline = time.monotonic() + max_total_wait_s
        counter = 0
        while True:
            counter += 1
            raw_path = Path(raw_status_dir) / f"status_{counter:05d}.txt" if raw_status_dir is not None else None
            result = self.check_status(session_number, raw_response_path=raw_path)
            if result.state in {"finished", "failed"}:
                return result
            if time.monotonic() >= deadline:
                raise SwissParamError(f"SwissParam session {session_number} exceeded max wait of {max_total_wait_s} s")
            self.sleeper(poll_interval_s)

    def retrieve(self, session_number: str, output_path: str | Path) -> Path:
        response = self._request("GET", f"{self.base_url}/retrievesession", params={"sessionNumber": session_number}, stream=True)
        target = Path(output_path)
        target.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = target.with_suffix(target.suffix + ".part")
        with tmp_path.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
        tmp_path.replace(target)
        if not tarfile.is_tarfile(target):
            raise SwissParamProtocolError(f"Retrieved payload is not a tar archive: {target.read_bytes()[:200]!r}")
        return target


def main() -> int:
    p = argparse.ArgumentParser(description="SwissParam command-line client")
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("health")
    submit = sub.add_parser("submit"); submit.add_argument("mol2", type=Path); submit.add_argument("--approach", default="both")
    status = sub.add_parser("status"); status.add_argument("session_number")
    retrieve = sub.add_parser("retrieve"); retrieve.add_argument("session_number"); retrieve.add_argument("output", type=Path)
    a = p.parse_args(); client = SwissParamClient()
    if a.command == "health": print(client.check_health())
    elif a.command == "submit": print(client.submit(a.mol2, approach=a.approach))
    elif a.command == "status": print(client.check_status(a.session_number).raw_text)
    else: print(client.retrieve(a.session_number, a.output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
