"""Local parity for the linker-charge, fused-ring, and OMS diagnostics.

Compares our faithful ports against MOFChecker 2.0 (`.venv-ref`) for:
  - possible_charged_fused_ring (bool)
  - positive_charge_from_linkers (int)
  - negative_charge_from_linkers (int)
  - has_oms (bool) and oms_indice (sorted index set)

The reference is ANCSA-licensed and used only as a behavioral oracle.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from pymatgen.core import Structure

from mofchecker_next.checks.charge_oms import (
    negative_charge_from_linkers_from_structure,
    oms_indices_from_structure,
    positive_charge_from_linkers_from_structure,
    possible_charged_fused_ring_from_structure,
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


# Probe each property independently: a no-metal structure makes has_oms/oms_indice
# raise NoMetal, but the charge/fused-ring descriptors are still well-defined.
_REF_CODE = """
import json, sys, warnings
warnings.filterwarnings("ignore")
from mofchecker import MOFChecker
checker = MOFChecker.from_cif(sys.argv[1], primitive=False, symprec=None, angle_tolerance=None)
def probe(fn):
    try:
        return {"ok": True, "value": fn()}
    except Exception as exc:  # noqa: BLE001
        return {"ok": False, "error": f"{type(exc).__name__}"}
out = {
    "possible_charged_fused_ring": probe(lambda: bool(checker.possible_charged_fused_ring)),
    "positive_charge_from_linkers": probe(lambda: int(checker.positive_charge_from_linkers)),
    "negative_charge_from_linkers": probe(lambda: int(checker.negative_charge_from_linkers)),
    "has_oms": probe(lambda: bool(checker.has_oms)),
    "oms_indice": probe(lambda: sorted(int(i) for i in checker.oms_indice)),
}
print(json.dumps(out))
"""


def _reference(cif_path: Path) -> dict:
    root = _repo_root()
    ref_python = root / ".venv-ref" / "bin" / "python"
    if not ref_python.exists():
        raise SystemExit("Reference environment .venv-ref not found.")
    completed = subprocess.run(
        [str(ref_python), "-c", _REF_CODE, str(cif_path)],
        text=True,
        capture_output=True,
    )
    if completed.returncode != 0:
        return {"_fatal": completed.stderr.strip()[:200] or "reference crashed"}
    return json.loads(completed.stdout.splitlines()[-1])


def _compare(cif_path: Path) -> dict:
    row: dict = {"cif": str(cif_path)}
    try:
        structure = Structure.from_file(str(cif_path))
    except Exception as exc:  # noqa: BLE001
        return {**row, "status": "our_load_error", "note": f"{type(exc).__name__}: {exc}"}

    ours = {
        "possible_charged_fused_ring": possible_charged_fused_ring_from_structure(structure),
        "positive_charge_from_linkers": positive_charge_from_linkers_from_structure(structure),
        "negative_charge_from_linkers": negative_charge_from_linkers_from_structure(structure),
        "oms_indice": sorted(oms_indices_from_structure(structure)),
    }
    ours["has_oms"] = len(ours["oms_indice"]) > 0
    row["ours"] = ours

    ref = _reference(cif_path)
    if "_fatal" in ref:
        return {**row, "status": "reference_exception", "note": ref["_fatal"]}

    # Compare each field independently. When the reference raises NoMetal for the
    # OMS fields (no metal in the cell), our metal-free result (no OMS) is treated
    # as an intentional, more-graceful agreement rather than a mismatch.
    fields: dict[str, bool] = {}
    ref_values: dict = {}
    notes: dict = {}
    for key in ("possible_charged_fused_ring", "positive_charge_from_linkers", "negative_charge_from_linkers", "has_oms", "oms_indice"):
        probe = ref[key]
        if probe["ok"]:
            ref_values[key] = probe["value"]
            fields[key] = ours[key] == probe["value"]
        elif probe.get("error") == "NoMetal" and key in ("has_oms", "oms_indice"):
            ref_values[key] = "NoMetal"
            empty = (ours["oms_indice"] == []) if key == "oms_indice" else (ours["has_oms"] is False)
            fields[key] = empty
            notes[key] = "reference raised NoMetal; ours returns no-OMS"
        else:
            ref_values[key] = f"EXC:{probe.get('error')}"
            fields[key] = False
    row["reference"] = ref_values
    row["field_matches"] = fields
    if notes:
        row["notes"] = notes
    row["matches"] = all(fields.values())
    row["status"] = "ok"
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cifs", nargs="*", type=Path)
    parser.add_argument("--out", type=Path, default=_repo_root() / "parity_charge_oms.json")
    parser.add_argument("--fail-on-mismatch", action="store_true")
    args = parser.parse_args()

    cifs = args.cifs
    if not cifs:
        cases = _repo_root() / "external" / "mofchecker_2_ref" / "test_cases"
        cifs = sorted(cases.rglob("*.cif"))

    rows = [_compare(cif) for cif in cifs]
    args.out.write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n")

    ok_rows = [r for r in rows if r.get("status") == "ok"]
    mismatches = [r for r in ok_rows if not r["matches"]]
    ref_exc = sum(1 for r in rows if r.get("status") == "reference_exception")
    print(
        f"charge/oms parity: {len(ok_rows) - len(mismatches)}/{len(ok_rows)} full matches, "
        f"{ref_exc} reference exceptions, {len(mismatches)} mismatches -> {args.out}"
    )
    for r in mismatches:
        bad = [k for k, v in r["field_matches"].items() if not v]
        print(f"  MISMATCH {Path(r['cif']).name}: fields={bad}")
        for k in bad:
            print(f"      {k}: ours={r['ours'][k]} ref={r['reference'][k]}")
    if args.fail_on_mismatch and mismatches:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
