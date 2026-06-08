"""Parity check on MODEL-GENERATED CIFs (not clean QMOFs), over the EVAL
geometric descriptor set on BOTH sides via the drop-in MOFChecker classes.

Generated CIFs across step dirs reuse basenames (sample_00000.cif, ...), so we
stage uniquely-named symlinks first to avoid id collisions.

  ~/projects/mofchecker-next/.venv/bin/python scripts/generated_parity.py \
      --dirs DIR1 DIR2 ... --per-dir 80 --workers 12 --out generated_parity_report.json
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qmof_parity import _repo_root  # noqa: E402

EVAL_GEO = [
    "has_carbon", "has_hydrogen", "has_metal",
    "has_atomic_overlaps",
    "has_overcoordinated_c", "has_overcoordinated_n", "has_overcoordinated_h",
    "has_undercoordinated_c", "has_undercoordinated_n",
    "has_lone_molecule", "has_suspicious_terminal_oxo",
    "has_geometrically_exposed_metal", "has_3d_connected_graph",
]


def _our_compute(path):
    out = {"id": os.path.basename(path)}
    try:
        from mofchecker_next import MOFChecker
        from pymatgen.core import Structure

        mc = MOFChecker(Structure.from_file(path))
    except Exception as exc:  # noqa: BLE001
        return {**out, "_load_error": f"{type(exc).__name__}"}
    for d in EVAL_GEO:
        try:
            out[d] = bool(getattr(mc, d))
        except Exception as exc:  # noqa: BLE001
            out[d] = {"__exc__": type(exc).__name__}
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dirs", nargs="+", required=True)
    ap.add_argument("--per-dir", type=int, default=0)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--out", type=Path, default=_repo_root() / "generated_parity_report.json")
    args = ap.parse_args()

    # Stage uniquely-named symlinks: <dirtag>__<basename>.cif
    stage = Path("/tmp/gen_parity_stage")
    if stage.exists():
        for p in stage.glob("*.cif"):
            p.unlink()
    stage.mkdir(parents=True, exist_ok=True)
    paths = []
    for d in args.dirs:
        dp = Path(d)
        tag = dp.name or dp.parent.name
        cifs = sorted(dp.glob("*.cif"))
        if args.per_dir > 0:
            cifs = cifs[: args.per_dir]
        for c in cifs:
            link = stage / f"{tag}__{c.name}"
            if not link.exists():
                link.symlink_to(c.resolve())
            paths.append(str(link))
    print(f"generated parity: {len(paths)} structures from {len(args.dirs)} dirs, {args.workers} workers")

    t0 = time.perf_counter()
    with Pool(args.workers) as pool:
        ours = {r["id"]: r for r in pool.map(_our_compute, paths)}
    print(f"  ours done in {time.perf_counter() - t0:.0f}s ({len(ours)} unique ids)")

    ref_python = _repo_root() / ".venv-ref" / "bin" / "python"
    worker = _repo_root() / "scripts" / "gen_ref_worker.py"
    paths_file = Path("/tmp/gen_parity_paths.json")
    ref_out = Path("/tmp/gen_parity_ref.jsonl")
    paths_file.write_text(json.dumps(paths))
    t0 = time.perf_counter()
    proc = subprocess.run([str(ref_python), str(worker), str(paths_file), str(ref_out), str(args.workers)], text=True)
    if proc.returncode != 0:
        raise SystemExit("reference worker failed")
    ref = {}
    for line in ref_out.read_text().splitlines():
        r = json.loads(line)
        ref[r["id"]] = r
    print(f"  reference done in {time.perf_counter() - t0:.0f}s ({len(ref)} unique ids)")

    stats = {d: {"match": 0, "compared": 0, "ref_exc": 0, "our_exc": 0} for d in EVAL_GEO}
    mismatches = []
    load_errors = 0
    for sid, o in ours.items():
        r = ref.get(sid)
        if r is None:
            continue
        if "_load_error" in o or "_load_error" in r:
            load_errors += 1
            continue
        for d in EVAL_GEO:
            ov, rv = o.get(d), r.get(d)
            if isinstance(rv, dict):
                stats[d]["ref_exc"] += 1; continue
            if isinstance(ov, dict):
                stats[d]["our_exc"] += 1; continue
            stats[d]["compared"] += 1
            if ov == rv:
                stats[d]["match"] += 1
            else:
                mismatches.append({"id": sid, "descriptor": d, "ours": ov, "ref": rv})

    report = {"n": len(paths), "load_errors": load_errors, "stats": stats, "mismatches": mismatches}
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    print(f"\n  load/compute errors skipped: {load_errors}")
    print(f"\n  {'descriptor':34s} {'match':>12s}  {'ref_exc':>7s} {'our_exc':>7s}")
    total_m = total_c = 0
    for d in EVAL_GEO:
        s = stats[d]
        total_m += s["match"]; total_c += s["compared"]
        flag = "" if s["match"] == s["compared"] else "  <-- MISMATCH"
        print(f"  {d:34s} {s['match']:>6d}/{s['compared']:<5d}  {s['ref_exc']:>7d} {s['our_exc']:>7d}{flag}")
    print(f"\n  OVERALL: {total_m}/{total_c} descriptor-comparisons match ({100*total_m/max(total_c,1):.3f}%)")
    print(f"  {len(mismatches)} mismatches -> {args.out}")
    for m in mismatches[:40]:
        print(f"    {m['id']:30s} {m['descriptor']:32s} ours={m['ours']} ref={m['ref']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
