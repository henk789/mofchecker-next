"""Benchmark mofchecker-next (Rust kernels) vs MOFChecker 2.0 (reference) on the
same generated CIFs and the same EVAL geometric descriptor set.

Runs new side in-process (this venv); old side via .venv-ref subprocess. Reports
wall time + throughput at 1 worker (pure per-structure speedup) and W workers.

  ~/projects/mofchecker-next/.venv/bin/python scripts/bench_speed.py \
      --dir <cif_dir> --n 120 --workers 12
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from qmof_parity import _repo_root  # noqa: E402

EVAL_GEO = [
    "has_carbon", "has_hydrogen", "has_metal", "has_atomic_overlaps",
    "has_overcoordinated_c", "has_overcoordinated_n", "has_overcoordinated_h",
    "has_undercoordinated_c", "has_undercoordinated_n", "has_lone_molecule",
    "has_suspicious_terminal_oxo", "has_geometrically_exposed_metal",
    "has_3d_connected_graph",
]


def time_new(paths, workers):
    from mofchecker_next.batch import check_structures
    t = time.perf_counter()
    res = check_structures(paths, n_workers=workers, descriptors=EVAL_GEO)
    return time.perf_counter() - t, sum(1 for r in res if "error" not in r)


def time_old(paths, workers):
    ref_python = _repo_root() / ".venv-ref" / "bin" / "python"
    worker = _repo_root() / "scripts" / "gen_ref_worker.py"
    pf = Path("/tmp/bench_paths.json"); pf.write_text(json.dumps(paths))
    of = Path("/tmp/bench_old.jsonl")
    t = time.perf_counter()
    p = subprocess.run([str(ref_python), str(worker), str(pf), str(of), str(workers)],
                       capture_output=True, text=True)
    dt = time.perf_counter() - t
    if p.returncode != 0:
        print(p.stderr[-500:]); raise SystemExit("old worker failed")
    ok = sum(1 for line in of.read_text().splitlines() if "_load_error" not in line)
    return dt, ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--n", type=int, default=120)
    ap.add_argument("--workers", type=int, default=12)
    args = ap.parse_args()

    paths = [str(p) for p in sorted(Path(args.dir).glob("*.cif"))[: args.n]]
    print(f"benchmark: {len(paths)} CIFs, descriptors={len(EVAL_GEO)}, workers={args.workers}\n")

    results = {}
    # 1-worker first (pure per-structure cost), then W workers.
    print("running new  (1 worker) ..."); results["new_1"] = time_new(paths, 1)
    print("running old  (1 worker) ..."); results["old_1"] = time_old(paths, 1)
    print(f"running new  ({args.workers} workers) ..."); results["new_w"] = time_new(paths, args.workers)
    print(f"running old  ({args.workers} workers) ..."); results["old_w"] = time_old(paths, args.workers)

    n = len(paths)
    def line(label, dt, ok):
        return f"  {label:24s} {dt:8.1f}s  {1000*dt/n:7.1f} ms/struct  {n/dt:7.1f} struct/s  (ok={ok})"

    print("\n=== RESULTS ===")
    print(line("OLD mofchecker (1w)", *results["old_1"]))
    print(line("NEW mofchecker-next (1w)", *results["new_1"]))
    print(line(f"OLD mofchecker ({args.workers}w)", *results["old_w"]))
    print(line(f"NEW mofchecker-next ({args.workers}w)", *results["new_w"]))
    s1 = results["old_1"][0] / results["new_1"][0]
    sw = results["old_w"][0] / results["new_w"][0]
    print(f"\n  speedup @1 worker:        {s1:5.2f}x")
    print(f"  speedup @{args.workers} workers:      {sw:5.2f}x")
    print(f"  new {args.workers}w vs old 1w (end-to-end): {results['old_1'][0]/results['new_w'][0]:5.2f}x")


if __name__ == "__main__":
    main()
