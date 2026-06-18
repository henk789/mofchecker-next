"""Reproducible speedup benchmark: mofchecker-next vs MOFChecker 2.0.

The structure set is *pinned* and persisted, so reruns (e.g. before/after your
optimizations) time the **same** structures -- otherwise timing comparisons are
meaningless. A fingerprint of the structure set is recorded; comparing against a
baseline that used a different set raises a loud warning.

Times three things per structure (median of --repeats, one warmup), on both
implementations, loading structures identically via pymatgen ``Structure.from_file``:
  - graph : build the "vesta" structure graph (same call both sides; isolates the
            pymatgen-version effect, which is NOT attributable to this project)
  - full  : construct the checker + read the diagnostic descriptor set (the
            end-to-end per-structure cost this project actually optimizes)
  - eqeq  : charge equilibration (ours = Rust kernel, ref = pyeqeq C++)

Usage:
  # first run -> selects + pins a structure set, writes a baseline report
  python scripts/benchmark_speedup.py --n 60 --out bench_before.json

  # after your changes -> reuses the SAME pinned structures, compares
  python scripts/benchmark_speedup.py --baseline bench_before.json --out bench_after.json

Set the CIF source with --source or $QMOF_DIR (defaults to the local QMOF
relaxed_structures checkout; falls back to the bundled reference test_cases).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import statistics
import subprocess
import sys
from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _default_source() -> Path:
    env = os.environ.get("QMOF_DIR")
    if env:
        return Path(env)
    qmof = Path("/projects/p2/p_fm_mofs/adit-mof-repro/data/qmof/raw/relaxed_structures")
    if qmof.is_dir():
        return qmof
    return _repo_root() / "external" / "mofchecker_2_ref" / "test_cases"


def _select_paths(source: Path, n: int, seed: int) -> list[str]:
    pool = sorted(str(p) for p in source.rglob("*.cif"))
    le150 = Path("/tmp/qmof_le150.json")
    if le150.exists():
        keep = set(json.load(open(le150)))
        filtered = [p for p in pool if Path(p).name in keep]
        if filtered:
            pool = filtered
    if not pool:
        raise SystemExit(f"No CIFs found under {source}")
    rng = random.Random(seed)
    rng.shuffle(pool)
    return pool[:n]


def _fingerprint(ours_rows: dict) -> str:
    """Hash the structure set (id + atom count) so two runs are comparable only
    if they used the identical structures."""
    items = sorted((sid, rec.get("n_atoms")) for sid, rec in ours_rows.items() if "n_atoms" in rec)
    return hashlib.sha256(repr(items).encode()).hexdigest()[:16]


def _run_worker(python: str, impl: str, paths: list[str], repeats: int, include_eqeq: bool) -> dict:
    worker = _repo_root() / "scripts" / "_benchmark_worker.py"
    paths_file = Path(f"/tmp/bench_paths_{impl}.json")
    out_file = Path(f"/tmp/bench_out_{impl}.json")
    paths_file.write_text(json.dumps(paths))
    cmd = [python, str(worker), impl, str(paths_file), str(out_file), str(repeats)]
    if include_eqeq:
        cmd.append("--eqeq")
    proc = subprocess.run(cmd, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"{impl} worker failed")
    return json.load(open(out_file))


def _aggregate(ours: dict, ref: dict, metrics: list[str]) -> dict:
    agg = {}
    common = [sid for sid in ours if sid in ref and "error" not in ours[sid] and "error" not in ref[sid]]
    for metric in metrics:
        o = [ours[sid][metric] for sid in common if metric in ours[sid] and metric in ref[sid]]
        r = [ref[sid][metric] for sid in common if metric in ours[sid] and metric in ref[sid]]
        if not o:
            continue
        per_struct_speedup = [rr / oo for oo, rr in zip(o, r) if oo > 0]
        agg[metric] = {
            "ours_ms_median": statistics.median(o),
            "ref_ms_median": statistics.median(r),
            "ours_ms_total": sum(o),
            "ref_ms_total": sum(r),
            "speedup_median": statistics.median(per_struct_speedup),
            "speedup_total": sum(r) / sum(o) if sum(o) else float("nan"),
            "n": len(o),
        }
    agg["_n_compared"] = len(common)
    return agg


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", type=Path, default=None, help="dir of CIFs (default: $QMOF_DIR or local QMOF)")
    parser.add_argument("--n", type=int, default=60)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--structures-file", type=Path, default=_repo_root() / "bench_structures.json",
                        help="pinned structure list; reused across runs unless --reselect")
    parser.add_argument("--reselect", action="store_true", help="re-pick structures even if the pin file exists")
    parser.add_argument("--no-eqeq", action="store_true", help="skip the EQeq timing")
    parser.add_argument("--baseline", type=Path, default=None, help="prior report to compare against")
    parser.add_argument("--out", type=Path, default=_repo_root() / "bench_report.json")
    args = parser.parse_args()

    # --- pinned structure set -------------------------------------------------
    if args.structures_file.exists() and not args.reselect:
        paths = json.loads(args.structures_file.read_text())
        print(f"Using pinned structure set from {args.structures_file} ({len(paths)} structures)")
    else:
        source = args.source or _default_source()
        paths = _select_paths(source, args.n, args.seed)
        args.structures_file.write_text(json.dumps(paths, indent=0))
        print(f"Selected {len(paths)} structures from {source} (seed={args.seed}) -> pinned to {args.structures_file}")

    include_eqeq = not args.no_eqeq
    metrics = ["graph", "full"] + (["eqeq"] if include_eqeq else [])

    print("Timing ours ...")
    ours = _run_worker(sys.executable, "ours", paths, args.repeats, include_eqeq)
    ref_python = _repo_root() / ".venv-ref" / "bin" / "python"
    if not ref_python.exists():
        raise SystemExit("Reference env .venv-ref not found.")
    print("Timing reference ...")
    ref = _run_worker(str(ref_python), "ref", paths, args.repeats, include_eqeq)

    fp = _fingerprint(ours)
    agg = _aggregate(ours, ref, metrics)
    report = {
        "fingerprint": fp,
        "n_structures": len(paths),
        "repeats": args.repeats,
        "aggregate": agg,
        "per_structure": {"ours": ours, "ref": ref},
    }
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    # --- report ---------------------------------------------------------------
    print(f"\nStructure-set fingerprint: {fp}   ({agg['_n_compared']} structures compared)")
    print(f"\n{'metric':8s} {'ours ms/str':>12s} {'ref ms/str':>12s} {'speedup(med)':>13s} {'speedup(tot)':>13s}")
    print("-" * 62)
    for m in metrics:
        if m not in agg:
            continue
        a = agg[m]
        print(f"{m:8s} {a['ours_ms_median']:12.2f} {a['ref_ms_median']:12.2f} "
              f"{a['speedup_median']:12.2f}x {a['speedup_total']:12.2f}x")
    print(f"\nNote: the `graph` row is dominated by the pymatgen version (same call both\n"
          f"sides) and is NOT attributable to this project; `full`/`eqeq` are.")
    print(f"-> {args.out}")

    # --- baseline comparison --------------------------------------------------
    if args.baseline:
        base = json.loads(Path(args.baseline).read_text())
        print(f"\n=== vs baseline {args.baseline} ===")
        if base.get("fingerprint") != fp:
            print(f"  !! WARNING: structure sets DIFFER (baseline {base.get('fingerprint')} vs now {fp}).")
            print(f"  !! Timing deltas are NOT comparable. Re-run the baseline with the same")
            print(f"  !! pinned --structures-file, or delete it and regenerate both runs.")
        else:
            print(f"  structure sets match ({fp}) -- deltas are valid.\n")
            bagg = base.get("aggregate", {})
            print(f"  {'metric':8s} {'ours now':>10s} {'ours base':>10s} {'delta':>9s} {'speedup now/base':>18s}")
            for m in metrics:
                if m in agg and m in bagg:
                    now = agg[m]["ours_ms_median"]; old = bagg[m]["ours_ms_median"]
                    d = (now - old) / old * 100 if old else float("nan")
                    print(f"  {m:8s} {now:10.2f} {old:10.2f} {d:+8.1f}% "
                          f"{agg[m]['speedup_median']:8.2f}x /{bagg[m]['speedup_median']:7.2f}x")
    return 0


if __name__ == "__main__":
    sys.exit(main())
