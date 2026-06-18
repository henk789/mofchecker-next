"""End-to-end batch throughput: mofchecker-next vs MOFChecker 2.0.

Times the wall-clock to run the geometric diagnostic set over a pinned set of
structures, at 1 and N worker processes, for both implementations (each loaded
via pymatgen Structure.from_file so only the implementation differs). Reports
ms/structure, structures/s, and the speedup -- the numbers for the README table.

Run on a dedicated node (e.g. via srun) for stable absolute numbers.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

# Full geometric diagnostic set (everything the default suite checks except the
# EQeq charge solve), including has_oms.
DESCRIPTORS_FULL = [
    "has_carbon", "has_hydrogen", "has_metal", "metal_number",
    "has_atomic_overlaps", "has_overcoordinated_c", "has_overcoordinated_n",
    "has_overcoordinated_h", "has_undercoordinated_c", "has_undercoordinated_n",
    "has_undercoordinated_rare_earth", "has_undercoordinated_alkali_alkaline",
    "has_lone_molecule", "has_3d_connected_graph", "has_suspicious_terminal_oxo",
    "has_geometrically_exposed_metal", "possible_charged_fused_ring",
    "positive_charge_from_linkers", "negative_charge_from_linkers", "has_oms",
]

# The accelerated subset used for generated-structure evaluation: the graph/
# geometry checks mofchecker-next actually speeds up. Excludes has_oms (identical
# pymatgen LocalStructOrderParams on both sides) and the cycle-charge checks.
DESCRIPTORS_EVAL = [
    "has_carbon", "has_hydrogen", "has_metal", "has_atomic_overlaps",
    "has_overcoordinated_c", "has_overcoordinated_n", "has_overcoordinated_h",
    "has_undercoordinated_c", "has_undercoordinated_n", "has_lone_molecule",
    "has_suspicious_terminal_oxo", "has_geometrically_exposed_metal",
    "has_3d_connected_graph",
]

DESCRIPTORS = DESCRIPTORS_FULL

_REF_CODE = """
import json, sys, time, warnings, os
from multiprocessing import Pool
warnings.filterwarnings("ignore")
paths = json.load(open(sys.argv[1])); workers = int(sys.argv[2]); desc = json.loads(sys.argv[3])
def _one(p):
    from pymatgen.core import Structure
    from mofchecker import MOFChecker
    try:
        mc = MOFChecker(Structure.from_file(p), None, symprec=None, angle_tolerance=None, primitive=False)
        for d in desc:
            try: getattr(mc, d)
            except Exception: pass
    except Exception: pass
# warmup one (imports, caches)
_one(paths[0])
t = time.perf_counter()
if workers == 1:
    for p in paths: _one(p)
else:
    with Pool(workers) as pool: pool.map(_one, paths)
print(time.perf_counter() - t)
"""


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _time_ours(paths, workers, descriptors) -> float:
    from mofchecker_next.batch import check_structures

    check_structures(paths[:1], n_workers=1, descriptors=descriptors)  # warmup
    t = time.perf_counter()
    check_structures(paths, n_workers=workers, descriptors=descriptors)
    return time.perf_counter() - t


def _time_ref(paths, workers, descriptors) -> float:
    ref_python = _repo_root() / ".venv-ref" / "bin" / "python"
    paths_file = Path("/tmp/throughput_paths.json")
    paths_file.write_text(json.dumps(paths))
    out = subprocess.run(
        [str(ref_python), "-c", _REF_CODE, str(paths_file), str(workers), json.dumps(descriptors)],
        text=True, capture_output=True,
    )
    if out.returncode != 0:
        raise SystemExit("ref throughput failed:\n" + out.stderr[-2000:])
    return float(out.stdout.strip().splitlines()[-1])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structures-file", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--set", choices=["full", "eval"], default="full",
                        help="'full' = all geometric diagnostics incl has_oms; "
                             "'eval' = accelerated subset (no has_oms / cycle-charge checks)")
    parser.add_argument("--out", type=Path, default=_repo_root() / "throughput_report.json")
    args = parser.parse_args()

    descriptors = DESCRIPTORS_EVAL if args.set == "eval" else DESCRIPTORS_FULL
    paths = json.loads(args.structures_file.read_text())
    n = len(paths)
    print(f"Throughput: {n} structures, descriptor set '{args.set}' ({len(descriptors)} descriptors), workers in [1, {args.workers}]\n")

    rows = {}
    for impl, fn in [("MOFChecker 2.0", _time_ref), ("mofchecker-next", _time_ours)]:
        rows[impl] = {}
        for w in (1, args.workers):
            dt = fn(paths, w, descriptors)
            rows[impl][w] = {"total_s": dt, "ms_per_struct": dt / n * 1e3, "per_s": n / dt}
            print(f"  {impl:18s} {w:2d} core(s): {dt:6.2f}s  {dt/n*1e3:7.1f} ms/struct  {n/dt:6.1f} struct/s")

    print(f"\n  {'config':28s} {'per structure':>14s} {'throughput':>12s} {'speedup':>9s}")
    for w in (1, args.workers):
        ref = rows["MOFChecker 2.0"][w]; ours = rows["mofchecker-next"][w]
        print(f"  {'MOFChecker 2.0 - ' + str(w) + ' core':28s} {ref['ms_per_struct']:11.1f} ms {ref['per_s']:9.1f} /s {'1x':>9s}")
        print(f"  {'mofchecker-next - ' + str(w) + ' core':28s} {ours['ms_per_struct']:11.1f} ms {ours['per_s']:9.1f} /s "
              f"{ref['ms_per_struct']/ours['ms_per_struct']:7.1f}x")

    args.out.write_text(json.dumps({"n_structures": n, "rows": rows}, indent=2) + "\n")
    print(f"\n-> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
