from __future__ import annotations

import argparse
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

TOOLS = {"Git": ["git"], "Conda": ["conda"], "Mamba": ["mamba"], "VMD": ["vmd"], "GROMACS": ["gmx", "gmx_mpi"], "ORCA": ["orca"]}


def _version(command: str) -> str:
    for flags in (["--version"], ["-version"], ["version"]):
        try:
            result = subprocess.run([command, *flags], capture_output=True, text=True, timeout=10, check=False)
        except (OSError, subprocess.TimeoutExpired):
            continue
        text = (result.stdout or result.stderr).strip()
        if text:
            return text.splitlines()[0][:300]
    return "version not resolved"


def build_inventory(scope_note: str) -> str:
    lines = ["# Environment inventory", "", f"- Generated UTC: `{datetime.now(timezone.utc).isoformat()}`", f"- Scope: `{scope_note}`", f"- OS: `{platform.platform()}`", f"- Python executable: `{sys.executable}`", f"- Python version: `{platform.python_version()}`", "", "## Tool discovery", "", "| Tool | Status | Executable | Version probe |", "|---|---|---|---|"]
    for label, commands in TOOLS.items():
        resolved = next(((command, shutil.which(command)) for command in commands if shutil.which(command)), None)
        if resolved is None:
            lines.append(f"| {label} | NOT FOUND | — | — |")
        else:
            command, executable = resolved
            lines.append(f"| {label} | FOUND | `{executable}` | `{_version(command)}` |")
    lines.extend(["", "> This is an environment inventory, not a request to reinstall working scientific software.", ""])
    return "\n".join(lines)


def main() -> int:
    p = argparse.ArgumentParser(description="Inventory local parameterization environment")
    p.add_argument("--output", type=Path, default=Path("reports/environment_inventory.md"))
    p.add_argument("--scope-note", default="machine where this command is executed")
    a = p.parse_args()
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(build_inventory(a.scope_note), encoding="utf-8")
    print(a.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
