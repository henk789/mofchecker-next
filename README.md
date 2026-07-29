<h1 align="center">
  mofchecker-next
</h1>

<p align="center">
  <b>A fast, drop-in replacement for <a href="https://github.com/lamalab-org/mofchecker">MOFChecker</a> 2.0</b> — same diagnostics, same API, <b>~12× faster</b> on the geometric diagnostics (~7× with open-metal-site detection), built on Rust kernels and rustworkx graph algorithms.
</p>

<p align="center">
    <a href="https://github.com/henk789/mofchecker-next/actions/workflows/release.yml">
        <img alt="CI" src="https://github.com/henk789/mofchecker-next/actions/workflows/release.yml/badge.svg" />
    </a>
    <a href="https://pypi.org/project/mofchecker-next/">
        <img alt="PyPI" src="https://img.shields.io/pypi/v/mofchecker-next" />
    </a>
    <a href="https://pypi.org/project/mofchecker-next/">
        <img alt="Python versions" src="https://img.shields.io/pypi/pyversions/mofchecker-next" />
    </a>
    <a href="https://github.com/henk789/mofchecker-next/blob/main/py/mofchecker_next/eqeq/LICENSE">
        <img alt="License" src="https://img.shields.io/pypi/l/mofchecker-next" />
    </a>
    <img alt="Speedup" src="https://img.shields.io/badge/vs%20MOFChecker%202.0-~12%C3%97%20faster-brightgreen" />
    <img alt="Built with Rust" src="https://img.shields.io/badge/built%20with-Rust%20%2B%20PyO3-orange" />
</p>

Designed for the workloads where the original is painful: **validating thousands of model-generated MOFs** (e.g. from a diffusion model), where the slow paths — floating-solvent extraction and dimensionality — dominate.

## 💪 Getting started

Install a published wheel with `pip install mofchecker-next`. Building from the
source distribution additionally requires Rust, a C compiler, and libclang
(`qhull-sys` generates bindings at build time); ordinary wheel installs do not.
Release CI installs and smoke-tests every built wheel before publishing.

```python
from mofchecker_next import MOFChecker

# Default: exact legacy MOFChecker 0.9.6 behavior used by published MOF metrics.
mc = MOFChecker.from_cif("structure.cif")     # also .from_ase(atoms) / MOFChecker(structure)
mc.has_atomic_overlaps, mc.has_lone_molecule, mc.has_oms, mc.metal_number
descriptors = mc.get_mof_descriptors()

# Newer behavior remains explicit:
modern = MOFChecker.from_cif("structure.cif", mode="2.0")
corrected = MOFChecker.from_cif("structure.cif", mode="next")
```

`mode="next"` currently resolves to provisional `next-dev-2`. Its named
`generated_desolvated_organic_mof_v1` composite requires input sanity, C/metal
scope, objective overlap/overcoordination/component checks, and true-element
Rust EQeq at declared total charge 0. Unvalidated C/N undercoordination,
hydrogen presence, exposed-metal/terminal-oxo/rare-earth rules, connectivity and
OMS remain report-only. The `|q| > 4.0` policy is still provisional, so this is a
shadow-production evaluator rather than frozen `next-1.0`. Results include the
exact composite fields, structured check evidence, source-tree/profile/input
manifests, and separate structure/local-chemistry/charge/component/scope
statuses. See `docs/NEXT_PRODUCTION_PLAN.md` for release gates and blockers.

`mode="0.9.6"` restores the legacy undercoordinated C/N rules, carbon-overcoordination
rule, rare-earth element set (no Sc/Y), `|q| > 3` EqEq threshold and EQeq's
CIF-label element parsing, unsplit floating-component flag, CIF parsing,
symmetry, and primitive defaults. It is the default; request `mode="2.0"`
explicitly for newer behavior.

Validating many structures? `mofchecker_next.batch` parallelizes across structures, builds each graph once, and never aborts the run on a single bad structure:

```python
from mofchecker_next.batch import check_structures

# inputs may be CIF paths, pymatgen Structures, ASE Atoms, or a mix
results = check_structures(inputs, n_workers=16)          # 0.9.6, all CPUs
modern_results = check_structures(inputs, n_workers=16, mode="2.0")
next_results = check_structures(inputs, n_workers=16, mode="next")
bad = [r for r in results if r["has_atomic_overlaps"]]

# subset to skip work: composition-only descriptors skip graph construction entirely
fast = check_structures(inputs, descriptors=["has_atomic_overlaps", "has_overcoordinated_c"])
```

The CLI also defaults to 0.9.6; use `mofchecker-next --mode next <CIFs...>` for corrected shadow scoring or `--mode 2.0` for compatibility auditing. Each result is a dict with `index`, `id`, `n_atoms`, and the requested descriptors. A structure that fails gets an `error` field (`on_error="record"`) instead of aborting the batch; an error in any composite field yields an indeterminate verdict and is invalid in `valid_rate_incl_errors`.

## 🚀 Installation

```bash
pip install mofchecker-next
```

Latest from source (needs a Rust toolchain):

```bash
pip install git+https://github.com/henk789/mofchecker-next.git
```

## ✨ Why use it

- ⚡ **Fast.** On 150-atom QMOF MOFs, the geometric diagnostic set runs in **~101 ms/structure** single-core vs **~1.20 s** for MOFChecker 2.0 (**~12×**); **~77 structures/s** across 10 cores. The win is the floating-solvent and 3D-connectivity paths, ported off networkx onto rustworkx; the numeric kernels (distances, contacts, connected components, OMS Voronoi/order-parameters, EQeq) are Rust.
- 🔌 **Drop-in.** `MOFChecker`-compatible class — same properties, same `get_mof_descriptors()`. Switch the import and existing code keeps working.
- ✅ **Legacy parity-verified.** 0 disagreements with MOFChecker 0.9.6 over all 325,984 full-QMOF descriptor comparisons. The 2.0 mode is near-parity, not yet a full exact port. See **Parity** below.
- 📦 **Built for batches.** Parallel `check_structures`, graph built once per structure, failures isolated.
- 🔋 **Bit-exact charges.** `has_high_charges` is a faithful Rust port of EQeq (bit-exact equilibrated charges).
- 🔁 **Reproducible.** `symmetry_hash` is deterministic (the reference's depends on Python hash randomization).

## ⚡ Performance

30 QMOF relaxed 150-atom MOFs, cold batch (each structure processed once), MOFChecker 2.0 vs mofchecker-next on identical inputs, dedicated compute node, 1 vs 10 cores. Reproduce with `scripts/benchmark_throughput.py` (the structure set is pinned + fingerprinted so before/after runs are comparable).

**Geometric diagnostic set** (the structural/graph checks; excludes open-metal-site detection):

| | per structure | throughput | speedup |
|---|---:|---:|---:|
| MOFChecker 2.0 — 1 core | 1201 ms | 0.8 /s | 1× |
| **mofchecker-next — 1 core** | **101 ms** | **9.9 /s** | **11.9×** |
| MOFChecker 2.0 — 10 cores | 162 ms | 6.2 /s | 1× |
| **mofchecker-next — 10 cores** | **13.1 ms** | **76.6 /s** | **12.4×** |

**Full geometric suite** (includes `has_oms` open-metal-site order parameters):

| | per structure | throughput | speedup |
|---|---:|---:|---:|
| MOFChecker 2.0 — 1 core | 1535 ms | 0.7 /s | 1× |
| **mofchecker-next — 1 core** | **140 ms** | **7.1 /s** | **10.9×** |
| MOFChecker 2.0 — 10 cores | 197 ms | 5.1 /s | 1× |
| **mofchecker-next — 10 cores** | **20.0 ms** | **49.9 /s** | **9.8×** |

Where the speedup comes from — it is concentrated in one check. MOFChecker's floating-solvent detection (`has_lone_molecule`) builds a 3×3×3 supercell graph via pymatgen `StructureGraph.__mul__` (networkx `union`/`relabel` of 27 copies, ~940 ms/structure here) and 3D-connectivity runs Larsen dimensionality over networkx. These are replaced by direct integer image-offset algorithms on a rustworkx graph — O(N+E), no supercell — making `has_lone_molecule` ~94× faster (~940 ms → ~10 ms); it dominates the reference's runtime, so it drives the speedup. `has_oms` now uses Rust for Voronoi facet-neighbor selection and pymatgen-compatible local order-parameter formulas. Speedup also grows with structure size and is higher on generated/distorted structures (more disconnected fragments).

## ⚙️ How it works

Python owns CIF/structure loading, pymatgen integration, and orchestration. The heavy lifting is delegated:

- **Rust** (`_rust` PyO3 extension): minimum-image distances, short contacts, neighbor candidates, connected components, graph degrees, OMS Voronoi/order-parameters, and the EQeq charge solve.
- **rustworkx** (`checks/_subgraph_rx.py`): floating-solvent / lone-molecule detection (finite connected components via an image-offset consistency test) and Larsen dimensionality (rank of the lattice-image vectors a component spans). These replace the networkx-heavy paths.
- **structuregraph_helpers** is retained for the logic-critical, non-hot pieces it does well: graph construction (tuned VESTA cutoffs) and the Weisfeiler–Lehman graph hashes.

The structure graph is built once per `MOFChecker` and reused across all checks.

## ✅ Parity

- **Legacy default (`mode="0.9.6"`):** 0 disagreements over 325,984 comparisons (all 20,374 QMOF CIFs × 16 shared Boolean diagnostics) against pinned `mofchecker==0.9.6`, with identical 16,259/20,374 = 79.8027% composite validity and zero errors.
- **2.0 mode:** exact on the original 250-real-QMOF and reference-CIF harnesses, but **not established as full parity**. The enriched 256-QMOF audit has 1/4096 descriptor disagreement (`has_lone_molecule`), an intentional corrected periodic-component result. A full 20,374-QMOF exact 2.0 audit has not been completed.
- **Generated (distorted) structures:** the same corrected floating-component algorithm can differ from 2.0 when a finite molecule wraps a cell boundary.

Reproduce: `scripts/qmof_parity.py` (real QMOFs), `scripts/generated_parity.py` (generated CIFs), `scripts/validate_subgraph_rx.py` (floating-solvent port). Point them at a local QMOF CIF directory with `QMOF_DIR=...`.

Parity says which reference a mode reproduces, not which reference is right. For
that, `scripts/undercoordination_review.py` builds a blinded, protocol-matched,
per-atom hand-labeling study of the C/N disagreements between 0.9.6 and 2.0 and
scores both versions against the labels — see
[`docs/UNDERCOORDINATION_ADJUDICATION.md`](docs/UNDERCOORDINATION_ADJUDICATION.md).
The requirements for a separately versioned corrected evaluator are listed in
[`docs/FIXED_CHECKER.md`](docs/FIXED_CHECKER.md); the staged implementation and
validation plan is [`docs/NEXT_MODE_PLAN.md`](docs/NEXT_MODE_PLAN.md).

## ⚠️ Limitations & deliberate differences

- **The default favors literature parity over speed.** `mode="0.9.6"` uses the
  reference's supercell floating-component algorithm; explicit `mode="2.0"`
  keeps the faster corrected implementation.
- **Healing not implemented.** `adding_hydrogen` / `adding_linker` raise `NotImplementedError`.
- **No porosity.** `is_porous` returns `None` (no bundled Zeo++).
- **`mode="2.0"` uses real elements for `has_high_charges`.** MOFChecker hands EQeq a
  pymatgen CIF and EQeq's parser reads two characters of the `<symbol><index>`
  label column, so one-letter elements (`C1`, `O2`, `U0`) miss its table and fall
  back to hydrogen's ionization row; structures where nothing matches score
  all-zero charges. `mode="2.0"` uses the correct elements; the default
  `mode="0.9.6"` reproduces the reference behaviour bit-exactly.
- **`has_lone_molecule` is *more* correct than the reference.** MOFChecker 2.0's supercell + in-cell-filter heuristic silently misses finite molecules that wrap the unit-cell boundary (the origin-cell copy is truncated at the supercell face). mofchecker-next detects them via a topological finite-component test. This is the only descriptor that ever disagrees with the reference, only on pathological/distorted structures (0 disagreements on real QMOFs).
- **Graph construction is still the floor.** pymatgen's VESTA neighbor-finding is unchanged; the speedup is in the graph *algorithms*, not bond perception.
- **Determinism.** `symmetry_hash` is deterministic by design and will not match the reference's randomized value across runs.

## ⚖️ License

The **published package is GPLv2**, because it bundles `py/mofchecker_next/eqeq/` — a faithful translation of [EQeq](https://github.com/lsmo-epfl/EQeq) (GPLv2, see `py/mofchecker_next/eqeq/LICENSE`) — and the GPL governs the combined work. The non-eqeq sources are **MIT** (`LICENSE`); for an MIT-only build, omit the `eqeq` subpackage and the `has_high_charges` diagnostic.

The MOFChecker 2.0 checkout used as the behavioral oracle (ANCSA 1.0) is **not** redistributed; see `external/REFERENCE.md` to reproduce it locally.

## 🛠️ For developers

<details>
<summary>Build, test, and release</summary>

```bash
python -m maturin develop --release          # build the Rust extension into the venv
python -m pytest -q                          # Python tests
cargo test --release --manifest-path rust/Cargo.toml   # Rust tests
# Source builds need libclang available for qhull-sys/bindgen.
```

**Layout**

- `py/mofchecker_next/` — Python package (`checks/`, `diagnostics.py`, the `eqeq` subpackage).
- `py/mofchecker_next/checks/_subgraph_rx.py` — rustworkx floating-solvent + dimensionality.
- `rust/` — the `_rust` PyO3 extension (geometry + OMS + EQeq kernels).
- `scripts/` — parity harnesses and the speed benchmark.
- `tests/` — Rust and Python unit tests.
- `docs/DIAGNOSTIC_INVENTORY.md` — per-diagnostic parity status.

**Making a release** — wheels are built and published by `.github/workflows/release.yml` via PyPI Trusted Publishing on a version tag:

```bash
git tag v0.1.0 && git push origin v0.1.0
```

</details>
