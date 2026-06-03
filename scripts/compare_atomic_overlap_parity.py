from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from pymatgen.core import Structure

from mofchecker_next.checks.geometry import (
    build_overlap_cutoff_matrix,
    check_atomic_overlaps,
    diagnostics_to_dicts,
    structure_to_arrays,
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _reference_json(cif_path: Path) -> dict:
    root = _repo_root()
    out = root / ".tmp-reference.json"
    subprocess.run(
        [sys.executable, str(root / "scripts" / "run_reference_mofchecker.py"), str(cif_path), "--out", str(out)],
        check=True,
    )
    payload = json.loads(out.read_text())
    out.unlink(missing_ok=True)
    return payload["reference"]


def _reference_radii(symbols: list[str]) -> dict[str, float]:
    root = _repo_root()
    ref_python = root / ".venv-ref" / "bin" / "python"
    code = """
import json
import sys
from mofchecker.checks.data import _get_covalent_radius

symbols = json.loads(sys.argv[1])
print(json.dumps({symbol: _get_covalent_radius(symbol) for symbol in symbols}, sort_keys=True))
"""
    completed = subprocess.run(
        [str(ref_python), "-c", code, json.dumps(sorted(set(symbols)))],
        check=True,
        text=True,
        capture_output=True,
    )
    return json.loads(completed.stdout)


def _cutoff_matrix_from_reference(structure: Structure) -> list[list[float]]:
    symbols = [site.specie.symbol for site in structure]
    radii = _reference_radii(symbols)
    _, atomic_numbers, _ = structure_to_arrays(structure)
    radii_by_z = {int(site.specie.Z): float(radii[site.specie.symbol]) for site in structure}
    return build_overlap_cutoff_matrix(atomic_numbers, radii_by_z)


def _involved_indices(diagnostics) -> list[int]:
    indices = set()
    for diagnostic in diagnostics:
        for atom in diagnostic.atoms:
            indices.add(atom.index)
    return sorted(indices)


def _expand_inputs(paths: list[Path], limit: int | None) -> list[Path]:
    root = _repo_root()
    expanded = []
    if not paths:
        paths = [root / "external" / "mofchecker_2_ref" / "test_cases"]
    for path in paths:
        if path.is_dir():
            expanded.extend(sorted(path.glob("**/*.cif")))
        else:
            expanded.append(path)
    if limit is not None:
        expanded = expanded[:limit]
    return expanded


def compare(cif_paths: list[Path]) -> list[dict]:
    rows = []
    for cif_path in cif_paths:
        try:
            reference = _reference_json(cif_path)
        except Exception as exc:
            rows.append({"cif": str(cif_path), "status": "reference_exception", "error": repr(exc), "matches": False})
            continue
        try:
            structure = Structure.from_file(cif_path)
            cutoff_matrix = _cutoff_matrix_from_reference(structure)
            diagnostics = check_atomic_overlaps(structure, cutoff_matrix)
        except Exception as exc:
            rows.append({"cif": str(cif_path), "status": "ours_exception", "error": repr(exc), "matches": False})
            continue
        reference_indices = sorted(int(index) for index in reference.get("overlapping_indices", []))
        our_indices = _involved_indices(diagnostics)
        reference_boolean = bool(reference["has_atomic_overlaps"])
        our_boolean = len(diagnostics) > 0
        row = {
                "cif": str(cif_path),
                "status": "ok",
                "reference_has_atomic_overlaps": reference_boolean,
                "ours_has_atomic_overlaps": our_boolean,
                "reference_overlapping_indices": reference_indices,
                "ours_overlapping_indices": our_indices,
                "boolean_matches": reference_boolean == our_boolean,
                "index_set_matches": reference_indices == our_indices,
                "contacts": len(diagnostics),
                "first_contacts": diagnostics_to_dicts(diagnostics[:5]),
            }
        row["matches"] = row["boolean_matches"] and row["index_set_matches"]
        rows.append(row)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare atomic-overlap diagnostics against MOFChecker 2.0.")
    parser.add_argument("paths", nargs="*", type=Path, help="CIF files or directories containing CIFs")
    parser.add_argument("--max", type=int, default=None, help="Maximum number of CIFs to compare")
    parser.add_argument("--json-out", type=Path, default=None, help="Write JSON report to this path")
    parser.add_argument("--fail-on-mismatch", action="store_true", help="Exit nonzero if any row mismatches or errors")
    args = parser.parse_args()

    rows = compare(_expand_inputs(args.paths, args.max))
    output = json.dumps(rows, indent=2, sort_keys=True)
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(output + "\n")
    print(output)
    if args.fail_on_mismatch and not all(row.get("matches", False) for row in rows):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
