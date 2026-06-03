"""Speed comparison: our reimplementation vs MOFChecker 2.0.

Runs each side in its own venv (ours in the current interpreter, the reference in
``.venv-ref``) and times comparable work with warmup + repeats, reporting the
median wall-clock per structure.

Three measurements per CIF:
  - graph build       : our structuregraph_helpers graph vs the reference's
  - eqeq charges      : our Rust kernel (end-to-end) vs pyeqeq's C++ (run_on_cif)
  - full diagnostics  : the full supported diagnostic set, graph built once
"""

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from pymatgen.core import Structure


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


# ---------------------------------------------------------------------------
# Our side (current venv).
# ---------------------------------------------------------------------------
def _our_timings(cif: str, repeats: int) -> dict:
    from mofchecker_next.checks import charge_oms as co
    from mofchecker_next.checks.graph import (
        build_structure_graph,
        floating_solvent_indices_from_structure,
        geometrically_exposed_metal_indices_from_graph,
    )
    from mofchecker_next.eqeq import compute_charges

    metals = co.METALS

    def time_it(fn, n):
        fn()  # warmup
        ts = []
        for _ in range(n):
            t0 = time.perf_counter()
            fn()
            ts.append(time.perf_counter() - t0)
        return statistics.median(ts)

    structure = Structure.from_file(cif)

    graph_t = time_it(lambda: build_structure_graph(structure), repeats)
    eqeq_t = time_it(lambda: compute_charges(structure), repeats)

    def full():
        graph = build_structure_graph(structure)
        co.oms_indices(structure, graph)
        co.positive_charge_indices(structure, graph)
        co.negative_charge_indices(structure, graph)
        co.fused_ring_indices(structure, graph)
        geometrically_exposed_metal_indices_from_graph(graph, metals)
        compute_charges(structure)

    full_t = time_it(full, repeats)
    return {"graph": graph_t, "eqeq": eqeq_t, "full": full_t, "n_atoms": len(structure)}


# ---------------------------------------------------------------------------
# Reference side (.venv-ref subprocess).
# ---------------------------------------------------------------------------
_REF_CODE = r"""
import json, sys, time, statistics, warnings, types, pathlib, tempfile, os
warnings.filterwarnings("ignore")
cif = sys.argv[1]; repeats = int(sys.argv[2])
# pyeqeq pkg_resources shim
base = pathlib.Path(sys.argv[3])
shim = types.ModuleType("pkg_resources")
shim.resource_filename = lambda p, n: str(base / n)
sys.modules["pkg_resources"] = shim
from mofchecker import MOFChecker
from pyeqeq.main import run_on_cif

def time_it(fn, n):
    fn()
    ts = []
    for _ in range(n):
        t0 = time.perf_counter(); fn(); ts.append(time.perf_counter() - t0)
    return statistics.median(ts)

def build_checker():
    c = MOFChecker.from_cif(cif, primitive=False, symprec=None, angle_tolerance=None)
    _ = c.graph
    return c

graph_t = time_it(build_checker, repeats)

# eqeq end-to-end as charge_check.py does it: structure -> cif -> run_on_cif
_c = build_checker()
def eqeq():
    with tempfile.NamedTemporaryFile("w", suffix=".cif", delete=False) as h:
        p = h.name
    try:
        _c.structure.to(filename=p, fmt="cif"); run_on_cif(p, verbose=False)
    finally:
        os.unlink(p)
eqeq_t = time_it(eqeq, repeats)

def full():
    c = MOFChecker.from_cif(cif, primitive=False, symprec=None, angle_tolerance=None)
    _ = c.graph
    try: c.has_oms
    except Exception: pass
    c.positive_charge_from_linkers
    c.negative_charge_from_linkers
    c.possible_charged_fused_ring
    c.has_geometrically_exposed_metal
    c.has_high_charges
full_t = time_it(full, repeats)
print(json.dumps({"graph": graph_t, "eqeq": eqeq_t, "full": full_t}))
"""


def _ref_timings(cif: str, repeats: int) -> dict:
    root = _repo_root()
    ref_python = root / ".venv-ref" / "bin" / "python"
    pyeqeq_dir = root / ".venv-ref" / "lib" / "python3.9" / "site-packages" / "pyeqeq"
    completed = subprocess.run(
        [str(ref_python), "-c", _REF_CODE, cif, str(repeats), str(pyeqeq_dir)],
        text=True,
        capture_output=True,
    )
    if completed.returncode != 0:
        return {"error": completed.stderr.strip()[-300:]}
    return json.loads(completed.stdout.splitlines()[-1])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cifs", nargs="*", type=str)
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args()

    cifs = args.cifs
    if not cifs:
        base = _repo_root() / "external" / "mofchecker_2_ref" / "test_cases"
        # Representative real periodic MOFs, ascending size.
        names = ["AFIPAH_clean.cif", "ALALUU_clean.cif", "GIMVAA_clean.cif", "FONQIJ_clean.cif"]
        cifs = [str(next(base.rglob(n))) for n in names]

    header = f"{'CIF':24s} {'atoms':>5s} | {'metric':10s} {'ours(ms)':>10s} {'ref(ms)':>10s} {'speedup':>8s}"
    print(header)
    print("-" * len(header))
    for cif in cifs:
        ours = _our_timings(cif, args.repeats)
        ref = _ref_timings(cif, args.repeats)
        name = Path(cif).name
        if "error" in ref:
            print(f"{name:24s} {ours['n_atoms']:5d} | reference error: {ref['error'][:60]}")
            continue
        for metric in ("graph", "eqeq", "full"):
            o = ours[metric] * 1e3
            r = ref[metric] * 1e3
            speed = r / o if o > 0 else float("nan")
            tag = name if metric == "graph" else ""
            atoms = f"{ours['n_atoms']:5d}" if metric == "graph" else "     "
            print(f"{tag:24s} {atoms} | {metric:10s} {o:10.2f} {r:10.2f} {speed:7.2f}x")
    return 0


if __name__ == "__main__":
    sys.exit(main())
