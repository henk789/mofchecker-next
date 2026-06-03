"""Reference-side QMOF descriptor batch worker (runs in .venv-ref).

Reads a JSON list of CIF paths on argv[1], computes the in-scope MOFChecker 2.0
descriptors for each (each probed independently so one failure does not mask the
rest), and writes JSONL results to argv[2]. Parallelized with multiprocessing.
"""

import json
import os
import sys
import warnings
from multiprocessing import Pool

warnings.filterwarnings("ignore")


def _compute(path):
    out = {"id": os.path.basename(path)}
    try:
        from mofchecker import MOFChecker
        from pymatgen.core import Structure

        # Load identically to our side (Structure.from_file) so the comparison
        # exercises the diagnostic LOGIC, not differences between CIF parsers
        # (MOFChecker.from_cif uses CifParser().get_structures(), which can yield
        # a subtly different structure than Structure.from_file and change
        # parser-sensitive checks like floating-solvent supercell extraction).
        structure = Structure.from_file(path)
        checker = MOFChecker(structure, None, symprec=None, angle_tolerance=None, primitive=False)
        _ = checker.graph
    except Exception as exc:  # noqa: BLE001
        return {**out, "_load_error": f"{type(exc).__name__}"}

    def probe(fn):
        try:
            return fn()
        except Exception as exc:  # noqa: BLE001
            return {"__exc__": type(exc).__name__}

    out["has_carbon"] = probe(lambda: bool(checker.has_carbon))
    out["has_hydrogen"] = probe(lambda: bool(checker.has_hydrogen))
    out["has_nitrogen"] = probe(lambda: bool(checker.checks["has_nitrogen"].is_ok))
    out["has_metal"] = probe(lambda: bool(checker.has_metal))
    out["metal_number"] = probe(lambda: int(checker.metal_number))
    out["has_overcoordinated_c"] = probe(lambda: bool(checker.has_overcoordinated_c))
    out["has_overcoordinated_n"] = probe(lambda: bool(checker.has_overcoordinated_n))
    out["has_undercoordinated_n"] = probe(lambda: bool(checker.has_undercoordinated_n))
    out["has_undercoordinated_rare_earth"] = probe(lambda: bool(checker.has_undercoordinated_rare_earth))
    out["has_undercoordinated_alkali_alkaline"] = probe(lambda: bool(checker.has_undercoordinated_alkali_alkaline))
    out["has_lone_molecule"] = probe(lambda: bool(checker.has_lone_molecule))
    out["has_3d_connected_graph"] = probe(lambda: bool(checker.has_3d_connected_graph))
    out["has_suspicious_terminal_oxo"] = probe(lambda: bool(checker.has_suspicious_terminal_oxo))
    out["has_geometrically_exposed_metal"] = probe(lambda: bool(checker.has_geometrically_exposed_metal))
    out["possible_charged_fused_ring"] = probe(lambda: bool(checker.possible_charged_fused_ring))
    out["positive_charge_from_linkers"] = probe(lambda: int(checker.positive_charge_from_linkers))
    out["negative_charge_from_linkers"] = probe(lambda: int(checker.negative_charge_from_linkers))
    out["has_oms"] = probe(lambda: bool(checker.has_oms))
    return out


def main():
    paths = json.load(open(sys.argv[1]))
    out_path = sys.argv[2]
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    with open(out_path, "w") as handle, Pool(workers) as pool:
        for i, row in enumerate(pool.imap_unordered(_compute, paths, chunksize=1)):
            handle.write(json.dumps(row) + "\n")
            handle.flush()
            if (i + 1) % 25 == 0:
                print(f"  reference: {i + 1}/{len(paths)}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
