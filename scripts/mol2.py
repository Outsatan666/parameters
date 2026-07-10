from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Iterable


class Mol2ParseError(ValueError):
    """Raised when a MOL2 file violates the minimum input contract."""


@dataclass