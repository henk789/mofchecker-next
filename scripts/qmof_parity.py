"""Large-scale parity check over the QMOF150 pool (QMOF structures <=150 atoms).

Runs our diagnostics (parallel, this venv) and the MOFChecker 2.0 reference
(parallel, .venv-ref) over a sample of QMOF CIFs and reports per-descriptor
parity. Our side is fast; the reference is the bottleneck, so both sides are
parallelized across structures.

The QMOF CIFs live in the sibling adit-mof-repro checkout (CC BY 4.0). Nothing is
copied into this repo.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import subprocess
import sys
import time
from multiprocessing import Pool
from pathlib import Path

QMOF_DIR = Path("/projects/p2/p_fm_mofs/adit-mof-repro/data/qmof/raw/relaxed_structures")

DESCRIPTORS = [
    "has_carbon", "has_hydrogen", "has_nitrogen", "has_metal", "metal_number",
    "has_overcoordinated_c", "has_overcoordinated_n", "has_undercoordinated_n",
    "has_undercoordinated_rare_earth", "has_undercoordinated_alkali_alkaline",
    "has_lone_molecule", "has_3d_connected_graph", "has_suspicious_terminal_oxo",
    "has_geometrically_exposed_metal", "possible_charged_fused_ring",
    "positive_charge_from_linkers", "negative_charge_from_linkers", "has_oms",
]


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _our_compute(path: str) -> dict:
    from pymatgen.core import Structure

    from mofchecker_next.checks import charge_oms as co
    from mofchecker_next.checks import composition as comp
    from mofchecker_next.checks import graph as g

    metals = co.METALS
    out = {"id": os.path.basename(path)}
    try:
        s = Structure.from_file(path)
    except Exception as exc:  # noqa: BLE001
        return {**out, "_load_error": f"{type(exc).__name__}"}
    try:
        out.update(comp.simple_global_diagnostics(s, metals))
        out["has_overcoordinated_c"] = len(g.overcoordinated_carbon_indices_from_structure(s, metals)) > 0
        out["has_overcoordinated_n"] = len(g.overcoordinated_nitrogen_indices_from_structure(s, metals)) > 0
        out["has_undercoordinated_n"] = len(g.undercoordinated_nitrogen_indices_from_structure(s, metals)) > 0
        out["has_undercoordinated_rare_earth"] = len(g.undercoordinated_rare_earth_indices_from_structure(s)) > 0
        out["has_undercoordinated_alkali_alkaline"] = len(g.undercoordinated_alkali_alkaline_indices_from_structure(s)) > 0
        out["has_lone_molecule"] = len(g.floating_solvent_indices_from_structure(s)) > 0
        out["has_3d_connected_graph"] = g.is_3d_connected_graph_from_structure(s)
        out["has_suspicious_terminal_oxo"] = len(g.false_oxo_indices_from_structure(s, metals)) > 0
        out["has_geometrically_exposed_metal"] = len(g.geometrically_exposed_metal_indices_from_structure(s, metals)) > 0
        out["possible_charged_fused_ring"] = co.possible_charged_fused_ring_from_structure(s)
        out["positive_charge_from_linkers"] = co.positive_charge_from_linkers_from_structure(s)
        out["negative_charge_from_linkers"] = co.negative_charge_from_linkers_from_structure(s)
        out["has_oms"] = co.has_oms_from_structure(s)
    except Exception as exc:  # noqa: BLE001
        return {**out, "_our_error": f"{type(exc).__name__}: {exc}"[:160]}
    return out


def _select_sample(n: int, seed: int) -> list[str]:
    pool = sorted(p.name for p in QMOF_DIR.glob("*.cif"))
    le150 = Path("/tmp/qmof_le150.json")
    if le150.exists():
        pool = sorted(json.load(open(le150)))
    rng = random.Random(seed)
    rng.shuffle(pool)
    return pool[:n]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=250)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--workers", type=int, default=min(20, (os.cpu_count() or 2) - 1))
    parser.add_argument("--out", type=Path, default=_repo_root() / "qmof_parity_report.json")
    args = parser.parse_args()

    sample = _select_sample(args.n, args.seed)
    paths = [str(QMOF_DIR / name) for name in sample]
    print(f"QMOF parity: {len(paths)} structures, {args.workers} workers")

    # Our side (parallel, this venv).
    t0 = time.perf_counter()
    with Pool(args.workers) as pool:
        ours = {r["id"]: r for r in pool.map(_our_compute, paths)}
    print(f"  ours done in {time.perf_counter() - t0:.0f}s")

    # Reference side (parallel, .venv-ref subprocess).
    ref_python = _repo_root() / ".venv-ref" / "bin" / "python"
    worker = _repo_root() / "scripts" / "_qmof_reference_worker.py"
    paths_file = Path("/tmp/qmof_parity_paths.json")
    ref_out = Path("/tmp/qmof_parity_ref.jsonl")
    paths_file.write_text(json.dumps(paths))
    t0 = time.perf_counter()
    proc = subprocess.run(
        [str(ref_python), str(worker), str(paths_file), str(ref_out), str(args.workers)],
        text=True,
    )
    if proc.returncode != 0:
        raise SystemExit("reference worker failed")
    ref = {}
    for line in ref_out.read_text().splitlines():
        r = json.loads(line)
        ref[r["id"]] = r
    print(f"  reference done in {time.perf_counter() - t0:.0f}s")

    # Compare.
    stats = {d: {"match": 0, "compared": 0, "ref_exc": 0} for d in DESCRIPTORS}
    mismatches = []
    load_errors = 0
    for sid, o in ours.items():
        r = ref.get(sid)
        if r is None:
            continue
        if "_load_error" in o or "_our_error" in o or "_load_error" in r:
            load_errors += 1
            continue
        for d in DESCRIPTORS:
            rv = r.get(d)
            ov = o.get(d)
            if isinstance(rv, dict) and "__exc__" in rv:
                # Reference raised. NoMetal on has_oms == our no-OMS (False) is a match.
                if rv["__exc__"] == "NoMetal" and d == "has_oms" and ov is False:
                    stats[d]["compared"] += 1
                    stats[d]["match"] += 1
                else:
                    stats[d]["ref_exc"] += 1
                continue
            stats[d]["compared"] += 1
            if ov == rv:
                stats[d]["match"] += 1
            else:
                mismatches.append({"id": sid, "descriptor": d, "ours": ov, "ref": rv})

    report = {"n": len(paths), "load_errors": load_errors, "stats": stats, "mismatches": mismatches}
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    print(f"\n  load/compute errors skipped: {load_errors}")
    print(f"\n  {'descriptor':36s} {'match':>12s}  {'ref_exc':>7s}")
    total_m = total_c = 0
    for d in DESCRIPTORS:
        s = stats[d]
        total_m += s["match"]; total_c += s["compared"]
        flag = "" if s["match"] == s["compared"] else "  <-- MISMATCH"
        print(f"  {d:36s} {s['match']:>6d}/{s['compared']:<5d}  {s['ref_exc']:>7d}{flag}")
    print(f"\n  OVERALL: {total_m}/{total_c} descriptor-comparisons match ({100*total_m/max(total_c,1):.3f}%)")
    print(f"  {len(mismatches)} mismatches -> {args.out}")
    for m in mismatches[:25]:
        print(f"    {m['id']:20s} {m['descriptor']:36s} ours={m['ours']} ref={m['ref']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
