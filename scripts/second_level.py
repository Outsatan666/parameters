from __future__ import annotations

import argparse
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

CANDIDATE_COMMANDS = ("ffparam", "FFParam", "ffparam-v2")


def discover_ffparam() -> tuple[str | None, str | None]:
    for command in CANDIDATE_COMMANDS:
        executable = shutil.which(command)
        if executable:
            result = subprocess.run([executable, "--help"], capture_output=True, text=True, timeout=30, check=False)
            return executable, (result.stdout or "") + (result.stderr or "")
    return None, None


def write_ffparam_installation_report(output: str | Path) -> Path:
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    executable, help_text = discover_ffparam()
    lines = ["# FFParam Tier-2 discovery report", "", f"- Generated UTC: `{datetime.now(timezone.utc).isoformat()}`", "- Preferred Tier-2 label from project specification: `FFParam-v2`", "- Scientific role: parameter optimization/validation, not automatic validation certification.", "", "## Discovery", ""]
    if executable is None:
        lines.extend(["- CLI availability: `NOT FOUND`", "- Version: `UNRESOLVED`", "- Source: `UNRESOLVED`", "- Installation method: `NOT PERFORMED`", "- Required QM backend: `UNRESOLVED`", "- Required MM backend: `UNRESOLVED`", "- License constraints: `UNRESOLVED`", "", "Tier 2 is intentionally blocked. No CLI flags are guessed from memory.", "Before enabling FFParam execution, inspect the authoritative source/manual and the actually installed package, then capture `<program> --help` output."])
    else:
        lines.extend(["- CLI availability: `FOUND`", f"- Executable: `{executable}`", "", "## Captured `--help`", "", "```text", help_text[:12000].rstrip(), "```", "", "Wrapper implementation must be based on this captured CLI."])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main() -> int:
    p = argparse.ArgumentParser(description="Tier-2 discovery; never guesses FFParam CLI")
    p.add_argument("--report", type=Path, default=Path("reports/ffparam_installation.md"))
    a = p.parse_args()
    print(write_ffparam_installation_report(a.report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
