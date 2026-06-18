"""Timing worker for benchmark_speedup.py. Runs under either venv.

Usage: python _benchmark_worker.py <impl> <paths.json> <out.json> <repeats> [--eqeq]
  impl = "ours"  -> times mofchecker_next
  impl = "ref"   -> times MOFChecker 2.0 (run under .venv-ref)

For each structure it times (median of <repeats>, one warmup):
  - graph     : build the structuregraph_helpers "vesta" graph (identical call
                on both sides; isolates the pymatgen-version effect)
  - full      : construct the checker and read the diagnostic descriptor set
                (graph built once internally) -- the end-to-end per-structure cost
  - eqeq      : charge-equilibration (ours: Rust kernel; ref: pyeqeq C++), optional

Both sides load the structure with pymatgen ``Structure.from_file`` so the only
difference is the implementation, not the CIF parser.
"""

import json
import statistics
import sys
import time
import warnings

warnings.filterwarnings("ignore")

# Diagnostic descriptors that both implementations expose as properties.
BENCH_DESCRIPTORS = [
    "has_carbon", "has_hydrogen", "has_metal", "metal_number",
    "has_atomic_overlaps", "has_overcoordinated_c", "has_overcoordinated_n",
    "has_overcoordinated_h", "has_undercoordinated_c", "has_undercoordinated_n",
    "has_undercoordinated_rare_earth", "has_undercoordinated_alkali_alkaline",
    "has_lone_molecule", "has_3d_connected_graph", "has_suspicious_terminal_oxo",
    "has_geometrically_exposed_metal", "possible_charged_fused_ring",
    "positive_charge_from_linkers", "negative_charge_from_linkers", "has_oms",
]


def _median_ms(fn, repeats):
    fn()  # warmup (also surfaces errors before timing)
    samples = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        samples.append(time.perf_counter() - t0)
    return statistics.median(samples) * 1e3


def _make_ours():
    from pymatgen.core import Structure
    from structuregraph_helpers.create import get_structure_graph

    from mofchecker_next import MOFChecker
    from mofchecker_next.eqeq import compute_charges

    def graph(s):
        return get_structure_graph(s, "vesta")

    def full(s):
        mc = MOFChecker(s)
        for name in BENCH_DESCRIPTORS:
            try:
                getattr(mc, name)
            except Exception:  # noqa: BLE001
                pass

    def eqeq(s):
        compute_charges(s)

    return Structure, graph, full, eqeq


def _make_ref():
    import os
    import pathlib
    import tempfile
    import types

    from pymatgen.core import Structure
    from structuregraph_helpers.create import get_structure_graph

    from mofchecker import MOFChecker

    def graph(s):
        return get_structure_graph(s, "vesta")

    def full(s):
        mc = MOFChecker(s, None, symprec=None, angle_tolerance=None, primitive=False)
        for name in BENCH_DESCRIPTORS:
            try:
                getattr(mc, name)
            except Exception:  # noqa: BLE001
                pass

    # pyeqeq with the removed-pkg_resources shim.
    base = pathlib.Path(__file__).resolve().parents[1] / ".venv-ref/lib/python3.9/site-packages/pyeqeq"
    shim = types.ModuleType("pkg_resources")
    shim.resource_filename = lambda pkg, name: str(base / name)
    sys.modules.setdefault("pkg_resources", shim)
    from pyeqeq.main import run_on_cif

    def eqeq(s):
        with tempfile.NamedTemporaryFile("w", suffix=".cif", delete=False) as handle:
            p = handle.name
        try:
            s.to(filename=p, fmt="cif")
            run_on_cif(p, verbose=False)
        finally:
            os.unlink(p)

    return Structure, graph, full, eqeq


def main():
    impl = sys.argv[1]
    paths = json.load(open(sys.argv[2]))
    out_path = sys.argv[3]
    repeats = int(sys.argv[4])
    do_eqeq = "--eqeq" in sys.argv[5:]

    Structure, graph, full, eqeq = _make_ours() if impl == "ours" else _make_ref()

    rows = {}
    for i, path in enumerate(paths):
        rec = {}
        try:
            structure = Structure.from_file(path)
            rec["n_atoms"] = len(structure)
            rec["graph"] = _median_ms(lambda: graph(structure), repeats)
            rec["full"] = _median_ms(lambda: full(structure), repeats)
            if do_eqeq:
                rec["eqeq"] = _median_ms(lambda: eqeq(structure), repeats)
        except Exception as exc:  # noqa: BLE001
            rec["error"] = f"{type(exc).__name__}: {exc}"[:160]
        rows[path.split("/")[-1]] = rec
        if (i + 1) % 10 == 0:
            print(f"  [{impl}] {i + 1}/{len(paths)}", file=sys.stderr, flush=True)
    json.dump(rows, open(out_path, "w"))


if __name__ == "__main__":
    main()
