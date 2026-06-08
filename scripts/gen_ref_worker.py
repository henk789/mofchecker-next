"""Reference-side worker (runs in .venv-ref / MOFChecker 2.0).

Reads a JSON list of CIF paths on argv[1], computes the EVAL geometric
descriptors for each (probed independently), writes JSONL to argv[2].
Loads via Structure.from_file to match our side's loader exactly.
"""
import json
import os
import sys
import warnings
from multiprocessing import Pool

warnings.filterwarnings("ignore")

EVAL_GEO = [
    "has_carbon", "has_hydrogen", "has_metal",
    "has_atomic_overlaps",
    "has_overcoordinated_c", "has_overcoordinated_n", "has_overcoordinated_h",
    "has_undercoordinated_c", "has_undercoordinated_n",
    "has_lone_molecule", "has_suspicious_terminal_oxo",
    "has_geometrically_exposed_metal", "has_3d_connected_graph",
]


def _compute(path):
    out = {"id": os.path.basename(path)}
    try:
        from mofchecker import MOFChecker
        from pymatgen.core import Structure

        structure = Structure.from_file(path)
        checker = MOFChecker(structure, None, symprec=None, angle_tolerance=None, primitive=False)
        _ = checker.graph
    except Exception as exc:  # noqa: BLE001
        return {**out, "_load_error": f"{type(exc).__name__}"}

    for d in EVAL_GEO:
        try:
            out[d] = bool(getattr(checker, d))
        except Exception as exc:  # noqa: BLE001
            out[d] = {"__exc__": type(exc).__name__}
    return out


def main():
    paths = json.load(open(sys.argv[1]))
    out_path = sys.argv[2]
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    with open(out_path, "w") as handle, Pool(workers) as pool:
        for i, row in enumerate(pool.imap_unordered(_compute, paths, chunksize=1)):
            handle.write(json.dumps(row) + "\n")
            handle.flush()
            if (i + 1) % 50 == 0:
                print(f"  reference: {i + 1}/{len(paths)}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
