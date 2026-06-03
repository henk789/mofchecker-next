from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


DESCRIPTORS = ("name", "formula", "has_atomic_overlaps", "has_metal")


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def run_reference(cif_path: Path, include_graph: bool = False) -> dict:
    root = _repo_root()
    ref_python = root / ".venv-ref" / "bin" / "python"
    if not ref_python.exists():
        raise SystemExit(
            "Reference environment not found. Create it with: "
            "uv venv --python 3.9 .venv-ref && "
            "uv pip install --python .venv-ref/bin/python -e external/mofchecker_2_ref numpy scipy pymatgen pandas"
        )

    code = """
import json
import sys
from mofchecker import MOFChecker

cif_path = sys.argv[1]
include_graph = sys.argv[2] == "1"
descriptors = ["name", "formula", "has_atomic_overlaps", "has_metal"]
checker = MOFChecker.from_cif(cif_path, primitive=False, symprec=None, angle_tolerance=None)
payload = checker.get_mof_descriptors(descriptors)
payload["overlapping_indices"] = checker.get_overlapping_indices()
payload["has_lone_molecule"] = checker.has_lone_molecule
payload["lone_molecule_indices"] = checker.lone_molecule_indices
payload["has_overcoordinated_c"] = checker.has_overcoordinated_c
payload["overcoordinated_c_indices"] = checker.overvalent_c_indices
payload["has_overcoordinated_n"] = checker.has_overcoordinated_n
payload["overcoordinated_n_indices"] = checker.checks["no_overcoordinated_nitrogen"].flagged_indices
payload["has_overcoordinated_h"] = checker.has_overcoordinated_h
payload["overcoordinated_h_indices"] = checker.overvalent_h_indices
payload["has_3d_connected_graph"] = checker.has_3d_connected_graph
payload["has_undercoordinated_c"] = checker.has_undercoordinated_c
payload["undercoordinated_c_indices"] = checker.undercoordinated_c_indices
payload["has_undercoordinated_n"] = checker.has_undercoordinated_n
payload["undercoordinated_n_indices"] = checker.undercoordinated_n_indices
payload["has_undercoordinated_rare_earth"] = checker.has_undercoordinated_rare_earth
payload["undercoordinated_rare_earth_indices"] = checker.undercoordinated_rare_earth_indices
payload["has_undercoordinated_alkali_alkaline"] = checker.has_undercoordinated_alkali_alkaline
payload["undercoordinated_alkali_alkaline_indices"] = checker.checks["no_undercoordinated_alkali_alkaline"].flagged_indices
payload["has_suspicious_terminal_oxo"] = checker.has_suspicious_terminal_oxo
payload["suspicious_terminal_oxo_indices"] = checker.suspicious_terminal_oxo_indices
payload["has_geometrically_exposed_metal"] = checker.has_geometrically_exposed_metal
payload["geometrically_exposed_metal_indices"] = checker.geometrically_exposed_metal_indice
from mofchecker.checks.data import _get_vdw_radius
payload["vdw_h_radius"] = float(_get_vdw_radius("H"))
payload["has_carbon"] = checker.has_carbon
payload["has_hydrogen"] = checker.has_hydrogen
payload["has_nitrogen"] = checker.checks["has_nitrogen"].is_ok
payload["metal_number"] = checker.metal_number

if include_graph:
    from mofchecker.definitions import METALS
    from mofchecker.checks.data import _get_covalent_radius
    graph = checker.graph.graph
    graph_edges = []
    edge_images = []
    for u, v, data in graph.edges(data=True):
        graph_edges.append([int(u), int(v)])
        edge_images.append([int(x) for x in data.get("to_jimage", (0, 0, 0))])
    payload["graph"] = {
        "atomic_numbers": [int(site.specie.Z) for site in checker.structure],
        "graph_edges": graph_edges,
        "edge_images": edge_images,
        "metal_atomic_numbers": sorted({int(element.Z) for element in checker.structure.composition.elements if element.symbol in METALS}),
        "metal_symbols": sorted(set(METALS)),
        "covalent_radii": {str(element.symbol): float(_get_covalent_radius(str(element.symbol))) for element in checker.structure.composition.elements},
    }
print(json.dumps(payload, sort_keys=True))
"""
    completed = subprocess.run(
        [str(ref_python), "-c", code, str(cif_path), "1" if include_graph else "0"],
        text=True,
        capture_output=True,
    )
    if completed.returncode != 0:
        raise SystemExit(completed.stderr or completed.stdout)
    return json.loads(completed.stdout)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run MOFChecker 2.0 reference descriptors.")
    parser.add_argument("cif", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--include-graph", action="store_true", help="Include reference-built graph arrays")
    args = parser.parse_args()

    payload = {
        "path": str(args.cif),
        "reference": run_reference(args.cif, include_graph=args.include_graph),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
