from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from pymatgen.core import Structure

from mofchecker_next.checks.composition import simple_global_diagnostics
from mofchecker_next.checks.geometry import overcoordinated_hydrogen_indices
from mofchecker_next.checks.graph import (
    floating_solvent_indices_from_structure,
    false_oxo_indices_from_structure,
    geometrically_exposed_metal_indices_from_structure,
    is_3d_connected_graph_from_structure,
    overcoordinated_carbon_indices_from_structure,
    overcoordinated_nitrogen_indices_from_structure,
    undercoordinated_carbon_indices_from_structure,
    undercoordinated_alkali_alkaline_indices_from_structure,
    undercoordinated_nitrogen_indices_from_structure,
    undercoordinated_rare_earth_indices_from_structure,
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _reference_json(cif_path: Path) -> dict:
    root = _repo_root()
    out = root / ".tmp-reference-graph.json"
    subprocess.run(
        [
            sys.executable,
            str(root / "scripts" / "run_reference_mofchecker.py"),
            str(cif_path),
            "--out",
            str(out),
            "--include-graph",
        ],
        check=True,
    )
    payload = json.loads(out.read_text())
    out.unlink(missing_ok=True)
    return payload["reference"]


def _flatten_indices(indices) -> list[int]:
    flattened = set()
    for item in indices:
        if isinstance(item, list):
            flattened.update(int(value) for value in item)
        else:
            flattened.add(int(item))
    return sorted(flattened)


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
            graph = reference["graph"]
            metal_symbols = graph["metal_symbols"]
            ours_global = simple_global_diagnostics(structure, graph["metal_symbols"])

            ours_c = overcoordinated_carbon_indices_from_structure(structure, metal_symbols)
            ours_n = overcoordinated_nitrogen_indices_from_structure(structure, metal_symbols)
            ours_lone_components = floating_solvent_indices_from_structure(structure)
            ours_3d = is_3d_connected_graph_from_structure(structure)
            ours_h = overcoordinated_hydrogen_indices(structure, reference["vdw_h_radius"])
            ours_under_c = undercoordinated_carbon_indices_from_structure(
                structure,
                metal_symbols,
                graph["covalent_radii"],
            )
            ours_under_n = undercoordinated_nitrogen_indices_from_structure(structure, metal_symbols)
            ours_under_re = undercoordinated_rare_earth_indices_from_structure(structure)
            ours_under_alk = undercoordinated_alkali_alkaline_indices_from_structure(structure)
            ours_false_oxo = false_oxo_indices_from_structure(structure, metal_symbols)
            ours_exposed_metal = geometrically_exposed_metal_indices_from_structure(structure, metal_symbols)
        except Exception as exc:
            rows.append({"cif": str(cif_path), "status": "ours_exception", "error": repr(exc), "matches": False})
            continue

        reference_lone_flat = _flatten_indices(reference["lone_molecule_indices"])
        ours_lone_flat = _flatten_indices(ours_lone_components)
        row = {
                "cif": str(cif_path),
                "status": "ok",
                "reference_has_lone_molecule": bool(reference["has_lone_molecule"]),
                "ours_has_lone_molecule_scaffold": len(ours_lone_components) > 0,
                "reference_lone_indices_flat": reference_lone_flat,
                "ours_lone_indices_flat": ours_lone_flat,
                "lone_boolean_matches": bool(reference["has_lone_molecule"]) == (len(ours_lone_components) > 0),
                "lone_index_set_matches": reference_lone_flat == ours_lone_flat,
                "reference_overcoordinated_c_indices": sorted(int(i) for i in reference["overcoordinated_c_indices"]),
                "ours_overcoordinated_c_indices": ours_c,
                "overcoordinated_c_matches": sorted(int(i) for i in reference["overcoordinated_c_indices"]) == ours_c,
                "reference_overcoordinated_n_indices": sorted(int(i) for i in reference["overcoordinated_n_indices"]),
                "ours_overcoordinated_n_indices": ours_n,
                "overcoordinated_n_matches": sorted(int(i) for i in reference["overcoordinated_n_indices"]) == ours_n,
                "reference_has_overcoordinated_h": bool(reference["has_overcoordinated_h"]),
                "reference_overcoordinated_h_indices": sorted(int(i) for i in reference["overcoordinated_h_indices"]),
                "ours_overcoordinated_h_indices": ours_h,
                "overcoordinated_h_matches": sorted(int(i) for i in reference["overcoordinated_h_indices"]) == ours_h,
                "reference_has_3d_connected_graph": bool(reference["has_3d_connected_graph"]),
                "ours_has_3d_connected_graph": ours_3d,
                "has_3d_connected_graph_matches": bool(reference["has_3d_connected_graph"]) == ours_3d,
                "reference_undercoordinated_c_indices": sorted(int(i) for i in reference["undercoordinated_c_indices"]),
                "ours_undercoordinated_c_indices": ours_under_c,
                "undercoordinated_c_matches": sorted(int(i) for i in reference["undercoordinated_c_indices"]) == ours_under_c,
                "reference_undercoordinated_n_indices": sorted(int(i) for i in reference["undercoordinated_n_indices"]),
                "ours_undercoordinated_n_indices": ours_under_n,
                "undercoordinated_n_matches": sorted(int(i) for i in reference["undercoordinated_n_indices"]) == ours_under_n,
                "reference_undercoordinated_rare_earth_indices": sorted(int(i) for i in reference["undercoordinated_rare_earth_indices"]),
                "ours_undercoordinated_rare_earth_indices": ours_under_re,
                "undercoordinated_rare_earth_matches": sorted(int(i) for i in reference["undercoordinated_rare_earth_indices"]) == ours_under_re,
                "reference_undercoordinated_alkali_alkaline_indices": sorted(int(i) for i in reference["undercoordinated_alkali_alkaline_indices"]),
                "ours_undercoordinated_alkali_alkaline_indices": ours_under_alk,
                "undercoordinated_alkali_alkaline_matches": sorted(int(i) for i in reference["undercoordinated_alkali_alkaline_indices"]) == ours_under_alk,
                "reference_has_suspicious_terminal_oxo": bool(reference["has_suspicious_terminal_oxo"]),
                "reference_suspicious_terminal_oxo_indices": sorted(int(i) for i in reference["suspicious_terminal_oxo_indices"]),
                "ours_suspicious_terminal_oxo_indices": ours_false_oxo,
                "suspicious_terminal_oxo_matches": sorted(int(i) for i in reference["suspicious_terminal_oxo_indices"]) == ours_false_oxo,
                "reference_has_geometrically_exposed_metal": bool(reference["has_geometrically_exposed_metal"]),
                "reference_geometrically_exposed_metal_indices": sorted(int(i) for i in reference["geometrically_exposed_metal_indices"]),
                "ours_geometrically_exposed_metal_indices": ours_exposed_metal,
                "geometrically_exposed_metal_matches": sorted(int(i) for i in reference["geometrically_exposed_metal_indices"]) == ours_exposed_metal,
                "simple_global_matches": {
                    key: reference[key] == ours_global[key]
                    for key in ("has_carbon", "has_hydrogen", "has_nitrogen", "has_metal", "metal_number")
                },
                "note": "C/N and lone molecule use structuregraph_helpers in our Python layer. H uses pymatgen radius-neighbor counting with the H VDW radius supplied by the reference wrapper.",
            }
        row["matches"] = (
            row["lone_boolean_matches"]
            and row["lone_index_set_matches"]
            and row["overcoordinated_c_matches"]
            and row["overcoordinated_n_matches"]
            and row["overcoordinated_h_matches"]
            and row["has_3d_connected_graph_matches"]
            and row["undercoordinated_c_matches"]
            and row["undercoordinated_n_matches"]
            and row["undercoordinated_rare_earth_matches"]
            and row["undercoordinated_alkali_alkaline_matches"]
            and row["suspicious_terminal_oxo_matches"]
            and row["geometrically_exposed_metal_matches"]
            and all(row["simple_global_matches"].values())
        )
        rows.append(row)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare graph-check scaffolding against MOFChecker 2.0 locally.")
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
