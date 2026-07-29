"""Immutable mode profiles.

Compatibility modes are frozen; ``next`` is a rolling alias for the explicitly
provisional corrected profile below.
"""

from __future__ import annotations

import hashlib
import importlib.util
import subprocess
from dataclasses import asdict, dataclass
from functools import cache
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from types import MappingProxyType


ADIT_PROBLEM_FLAGS = (
    "has_atomic_overlaps", "has_overcoordinated_c", "has_overcoordinated_n",
    "has_overcoordinated_h", "has_undercoordinated_c", "has_undercoordinated_n",
    "has_undercoordinated_rare_earth", "has_undercoordinated_alkali_alkaline",
    "has_stray_atom", "has_lone_molecule", "has_suspicious_terminal_oxo",
    "has_high_charges", "has_geometrically_exposed_metal",
)
ADIT_PRESENCE_FLAGS = ("has_carbon", "has_hydrogen", "has_metal")
LEGACY_PROBLEM_FLAGS = tuple(
    flag for flag in ADIT_PROBLEM_FLAGS if flag != "has_stray_atom"
)
NEXT_COMPOSITE_NAME = "generated_desolvated_organic_mof_v1"
NEXT_PROBLEM_FLAGS = (
    "has_atomic_overlaps", "has_overcoordinated_c", "has_overcoordinated_n",
    "has_overcoordinated_h", "has_stray_atom", "has_lone_molecule",
    "has_high_charges",
)
NEXT_PRESENCE_FLAGS = ("input_sanity", "has_carbon", "has_metal")
NEXT_REPORT_ONLY_FIELDS = (
    "has_hydrogen", "has_undercoordinated_c", "has_undercoordinated_n",
    "has_undercoordinated_rare_earth", "has_undercoordinated_alkali_alkaline",
    "has_3d_connected_graph", "has_suspicious_terminal_oxo",
    "has_geometrically_exposed_metal", "has_oms", "possible_charged_fused_ring",
    "positive_charge_from_linkers", "negative_charge_from_linkers",
)
NEXT_STATUS_GROUPS = (
    ("structure_status", "problem", ("has_atomic_overlaps", "has_overcoordinated_h")),
    ("local_chemistry_status", "problem", (
        "has_overcoordinated_c", "has_overcoordinated_n",
    )),
    ("charge_status", "problem", ("has_high_charges",)),
    ("component_status", "problem", ("has_stray_atom", "has_lone_molecule")),
    ("scope_status", "presence", ("has_carbon", "has_metal")),
)


@dataclass(frozen=True)
class ModeProfile:
    id: str
    descriptor_set: str
    composite_name: str = "adit_mofasa_all_passed"
    problem_flags: tuple[str, ...] = ADIT_PROBLEM_FLAGS
    presence_flags: tuple[str, ...] = ADIT_PRESENCE_FLAGS
    status_groups: tuple[tuple[str, str, tuple[str, ...]], ...] = ()
    report_only_fields: tuple[str, ...] = ()
    legacy_input: bool = False
    corrected_carbon_images: bool = False
    corrected_input: bool = False
    total_charge: float | None = None
    eqeq_threshold: float = 4.0
    eqeq_charge_sum_tolerance: float = 1e-9
    provisional: bool = False


_PROFILES = {
    "0.9.6": ModeProfile(
        "0.9.6", "legacy", problem_flags=LEGACY_PROBLEM_FLAGS,
        legacy_input=True, eqeq_threshold=3.0,
    ),
    "2.0": ModeProfile("2.0", "modern"),
    "next": ModeProfile(
        "next-dev-2",
        "next",
        composite_name=NEXT_COMPOSITE_NAME,
        problem_flags=NEXT_PROBLEM_FLAGS,
        presence_flags=NEXT_PRESENCE_FLAGS,
        status_groups=NEXT_STATUS_GROUPS,
        report_only_fields=NEXT_REPORT_ONLY_FIELDS,
        corrected_carbon_images=True,
        corrected_input=True,
        total_charge=0.0,
        provisional=True,
    ),
}
PROFILES = MappingProxyType(_PROFILES)
MODES = tuple(PROFILES)
DEFAULT_MODE = "0.9.6"


def resolve_profile(mode: str) -> ModeProfile:
    try:
        return PROFILES[mode]
    except KeyError:
        raise ValueError(f"mode must be one of {MODES}, got {mode!r}") from None


@cache
def source_identity() -> dict:
    """Version and content identity; the tree digest also identifies dirty builds."""
    package_dir = Path(__file__).resolve().parent
    root = package_dir.parents[1] if package_dir.parent.name == "py" else package_dir.parent
    files = sorted(package_dir.rglob("*.py"))
    files += sorted(path for path in package_dir.glob("_rust*") if path.is_file())
    files += sorted((root / "rust").rglob("*.rs"))
    files += [path for path in (root / "pyproject.toml", root / "rust" / "Cargo.toml") if path.exists()]
    digest = hashlib.sha256()
    for path in files:
        try:
            identity_path = path.relative_to(root).as_posix()
        except ValueError:
            identity_path = path.name
        digest.update(identity_path.encode() + b"\0")
        digest.update(path.read_bytes())
    rust_spec = importlib.util.find_spec("mofchecker_next._rust")
    rust_path = Path(rust_spec.origin) if rust_spec and rust_spec.origin else None
    rust_extension_sha256 = (
        hashlib.sha256(rust_path.read_bytes()).hexdigest()
        if rust_path and rust_path.is_file() else None
    )
    try:
        package_version = version("mofchecker-next")
    except PackageNotFoundError:
        package_version = "unknown"
    try:
        if package_dir.parent.name != "py" or not (root / ".git").exists():
            raise OSError("not an editable source checkout")
        source_sha = subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"], text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
        source_dirty = bool(subprocess.check_output(
            ["git", "-C", str(root), "status", "--porcelain"], text=True,
            stderr=subprocess.DEVNULL,
        ).strip())
    except (OSError, subprocess.CalledProcessError):
        source_sha, source_dirty = None, None
    return {
        "package_version": package_version,
        "source_sha": source_sha,
        "source_dirty": source_dirty,
        "source_tree_sha256": digest.hexdigest(),
        "rust_extension_sha256": rust_extension_sha256,
    }


def profile_manifest(mode: str) -> dict:
    profile = resolve_profile(mode)
    return {
        **source_identity(),
        "requested_mode": mode,
        "resolved_profile": profile.id,
        "profile": asdict(profile),
    }
