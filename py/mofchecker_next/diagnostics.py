from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal


@dataclass(frozen=True)
class AtomRef:
    index: int
    image: tuple[int, int, int] = (0, 0, 0)


@dataclass(frozen=True)
class Diagnostic:
    check: str
    severity: Literal["info", "warning", "error"]
    message: str
    atoms: list[AtomRef]
    values: dict = field(default_factory=dict)
    suggested_action: str | None = None
