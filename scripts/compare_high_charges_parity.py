"""Local parity for the EQeq ``has_high_charges`` diagnostic.

Compares our Rust EQeq port against the EQeq oracle (``pyeqeq``) bundled in
``.venv-ref``. Our kernel is a faithful translation of EQeq's C++ and matches it
bit-for-bit at the 3-decimal output precision.

Important: EQeq ships a hand-rolled CIF parser that reads the element symbol from
the *label* column and keeps two characters. Two-letter elements survive
(``Pr1`` -> ``Pr``); single-letter labels (``C1``, ``O1``, ``H1``) are not in its
table, so they silently fall back to hydrogen. MOFChecker 2.0's ``has_high_charges``
feeds EQeq a pymatgen-written CIF whose label column is ``C1``/``O1``/..., so the
reference computes charges treating every single-letter element as hydrogen -- a
reference-side parsing bug, not an algorithm difference.

To measure true algorithm parity we therefore hand the oracle an "EQeq-safe" CIF
whose label column equals the element symbol, so EQeq parses the same structure
our kernel does. Both the per-atom charges and the diagnostic boolean are checked.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from pymatgen.core import Structure

from mofchecker_next.eqeq import (
    DEFAULT_CHARGE_THRESHOLD,
    compute_charges,
    high_charge_indices,
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


# Runs in .venv-ref: shims the removed pkg_resources, then calls pyeqeq. The CIF
# passed in has already been written by pymatgen, exactly mirroring MOFChecker
# 2.0's charge_check.py (structure.to(cif) -> run_on_cif), so both sides consume
# the same structure and only EQeq's own parse/algorithm is exercised.
_ORACLE_CODE = """
import json, sys, types, pathlib, warnings
warnings.filterwarnings("ignore")
base = pathlib.Path(sys.argv[2])
shim = types.ModuleType("pkg_resources")
shim.resource_filename = lambda pkg, name: str(base / name)
sys.modules["pkg_resources"] = shim
from pyeqeq.main import run_on_cif
try:
    charges = run_on_cif(sys.argv[1], verbose=False)
    print(json.dumps({"ok": True, "charges": list(charges)}))
except Exception as exc:  # noqa: BLE001 - report any oracle failure as data
    print(json.dumps({"ok": False, "error": f"{type(exc).__name__}: {exc}"}))
"""


def _eqeq_safe_cif(structure) -> str:
    """A minimal CIF whose label column equals the element symbol, so EQeq's
    crude parser reads the correct element for single-letter species too."""
    a, b, c, alpha, beta, gamma = structure.lattice.parameters
    lines = [
        "data_eqeq",
        f"_cell_length_a {a:.10f}",
        f"_cell_length_b {b:.10f}",
        f"_cell_length_c {c:.10f}",
        f"_cell_angle_alpha {alpha:.10f}",
        f"_cell_angle_beta {beta:.10f}",
        f"_cell_angle_gamma {gamma:.10f}",
        "loop_",
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
    ]
    for site in structure:
        sym = site.specie.symbol
        f = site.frac_coords
        lines.append(f"{sym} {sym} {f[0]:.10f} {f[1]:.10f} {f[2]:.10f}")
    return "\n".join(lines) + "\n"


def _oracle_charges(structure) -> dict:
    import tempfile

    root = _repo_root()
    ref_python = root / ".venv-ref" / "bin" / "python"
    pyeqeq_dir = root / ".venv-ref" / "lib" / "python3.9" / "site-packages" / "pyeqeq"
    if not ref_python.exists():
        raise SystemExit("Reference environment .venv-ref not found.")
    with tempfile.NamedTemporaryFile("w", suffix=".cif", delete=False) as handle:
        handle.write(_eqeq_safe_cif(structure))
        cif_path = Path(handle.name)
    try:
        completed = subprocess.run(
            [str(ref_python), "-c", _ORACLE_CODE, str(cif_path), str(pyeqeq_dir)],
            text=True,
            capture_output=True,
        )
    finally:
        cif_path.unlink(missing_ok=True)
    if completed.returncode != 0:
        return {"ok": False, "error": completed.stderr.strip() or "oracle crashed"}
    return json.loads(completed.stdout)


def _compare_cif(cif_path: Path, threshold: float) -> dict:
    row: dict = {"cif": str(cif_path)}
    try:
        structure = Structure.from_file(str(cif_path))
    except Exception as exc:  # noqa: BLE001
        row["status"] = "our_load_error"
        row["note"] = f"{type(exc).__name__}: {exc}"
        return row

    our_charges = compute_charges(structure)
    our_indices = high_charge_indices(structure, threshold=threshold)
    row["ours_has_high_charges"] = bool(our_indices)
    row["ours_high_charge_indices"] = our_indices
    row["ours_max_abs_charge"] = float(np.abs(our_charges).max()) if len(our_charges) else 0.0

    oracle = _oracle_charges(structure)
    if not oracle["ok"]:
        row["status"] = "reference_exception"
        row["note"] = oracle["error"][:200]
        return row

    ref_charges = np.asarray(oracle["charges"], dtype=float)
    ref_high = bool((np.abs(ref_charges) > threshold).any())
    row["reference_has_high_charges"] = ref_high
    row["reference_max_abs_charge"] = float(np.abs(ref_charges).max()) if len(ref_charges) else 0.0
    row["boolean_matches"] = row["ours_has_high_charges"] == ref_high
    if len(ref_charges) != len(our_charges):
        row["status"] = "atom_count_mismatch"
        row["note"] = f"ours {len(our_charges)} vs oracle {len(ref_charges)}"
        return row
    # Per-atom charge agreement at EQeq's 3-decimal output precision.
    max_abs_dev = float(np.abs(our_charges - ref_charges).max()) if len(our_charges) else 0.0
    row["max_abs_charge_deviation"] = max_abs_dev
    row["charges_match"] = max_abs_dev < 1e-3
    row["status"] = "ok"
    row["matches"] = row["boolean_matches"] and row["charges_match"]
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "cifs",
        nargs="*",
        type=Path,
        help="CIF files (default: all reference test_cases CIFs).",
    )
    parser.add_argument("--out", type=Path, default=_repo_root() / "parity_high_charges.json")
    parser.add_argument("--threshold", type=float, default=DEFAULT_CHARGE_THRESHOLD)
    parser.add_argument("--fail-on-mismatch", action="store_true")
    args = parser.parse_args()

    cifs = args.cifs
    if not cifs:
        cases = _repo_root() / "external" / "mofchecker_2_ref" / "test_cases"
        cifs = sorted(cases.rglob("*.cif"))

    rows = [_compare_cif(cif, args.threshold) for cif in cifs]
    args.out.write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n")

    ok_rows = [r for r in rows if r.get("status") == "ok"]
    mismatches = [r for r in ok_rows if not r.get("matches")]
    worst = max((r.get("max_abs_charge_deviation", 0.0) for r in ok_rows), default=0.0)
    exact = sum(1 for r in ok_rows if r.get("charges_match"))
    print(
        f"high_charges parity: {len(ok_rows)} compared, {exact}/{len(ok_rows)} "
        f"bit-exact charges, worst max|Δq|={worst:.6f} -> {args.out}"
    )
    for r in mismatches:
        print(
            f"  MISMATCH {Path(r['cif']).name}: "
            f"bool ours={r.get('ours_has_high_charges')} ref={r.get('reference_has_high_charges')}, "
            f"max|Δq|={r.get('max_abs_charge_deviation')}"
        )
    if args.fail_on_mismatch and mismatches:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
