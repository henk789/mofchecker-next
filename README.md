# mofchecker-next

**A fast, drop-in replacement for [MOFChecker](https://github.com/lamalab-org/mofchecker) 2.0** — same diagnostics, same API, **~20–25× faster per core**, built on Rust kernels and rustworkx graph algorithms.

Designed for the workloads where the original is painful: **validating thousands of model-generated MOFs** (e.g. from a diffusion model), where the slow paths — floating-solvent extraction and dimensionality — dominate.

```python
from mofchecker_next import MOFChecker

mc = MOFChecker.from_cif("structure.cif")
mc.has_atomic_overlaps, mc.has_lone_molecule, mc.has_oms, mc.metal_number
descriptors = mc.get_mof_descriptors()        # OrderedDict of every diagnostic
```

---

## ✨ Why use it

- ⚡ **Fast.** On 150-atom MOFs, **57.8 ms/structure** single-core vs **1.34 s** for MOFChecker 2.0 (**23×**); **149 structures/s** across 10 cores. The hot paths (floating solvent, 3D-connectivity) were ported off networkx; the numeric kernels (distances, contacts, connected components, EQeq) are Rust.
- 🔌 **Drop-in.** `MOFChecker`-compatible class — same properties, same `get_mof_descriptors()`. Switch the import and existing code keeps working.
- ✅ **Parity-verified.** 100% agreement with MOFChecker 2.0 on real QMOFs (4500/4500 descriptor-comparisons over 250 structures × 18 descriptors; 16/16 on the reference test CIFs). See **Parity** below.
- 📦 **Built for batches.** `mofchecker_next.batch` parallelizes across structures, builds each graph once, and never aborts the run on a single bad structure.
- 🔋 **Bit-exact charges.** `has_high_charges` is a faithful Rust port of EQeq (bit-exact equilibrated charges).
- 🔁 **Reproducible.** `symmetry_hash` is deterministic (the reference's depends on Python hash randomization).

## ⚡ Performance

120 generated 150-atom CIFs, full geometric descriptor set, MOFChecker 2.0 vs mofchecker-next on identical inputs (10-core node):

| | per structure | throughput | speedup |
|---|---:|---:|---:|
| MOFChecker 2.0 — 1 core | 1338 ms | 0.7 /s | 1× |
| **mofchecker-next — 1 core** | **57.8 ms** | **17.3 /s** | **23.2×** |
| MOFChecker 2.0 — 10 cores | 174.9 ms | 5.7 /s | — |
| **mofchecker-next — 10 cores** | **6.7 ms** | **149.3 /s** | **26.1×** |

Where the speedup comes from: MOFChecker's floating-solvent check builds a 3×3×3 supercell graph via pymatgen `StructureGraph.__mul__` (networkx `union`/`relabel` of 27 copies, ~2.3 s/structure), and 3D-connectivity runs Larsen dimensionality over networkx. Both are replaced by direct integer image-offset algorithms on a rustworkx graph — O(N+E), no supercell — while the geometry kernels run in Rust.

## ⚙️ How it works

Python owns CIF/structure loading, pymatgen integration, and orchestration. The heavy lifting is delegated:

- **Rust** (`_rust` PyO3 extension): minimum-image distances, short contacts, neighbor candidates, connected components, graph degrees, and the EQeq charge solve.
- **rustworkx** (`checks/_subgraph_rx.py`): floating-solvent / lone-molecule detection (finite connected components via an image-offset consistency test) and Larsen dimensionality (rank of the lattice-image vectors a component spans). These replace the networkx-heavy paths in `structuregraph_helpers` / pymatgen.
- **structuregraph_helpers** is still used for the logic-critical, non-hot pieces it does well: graph construction (tuned VESTA cutoffs) and the Weisfeiler–Lehman graph hashes.

The structure graph is built once per `MOFChecker` and reused across all checks.

## 📦 Batch validation

```python
from mofchecker_next.batch import check_structures, check_structure

# inputs may be CIF paths, pymatgen Structures, ASE Atoms, or a mix
results = check_structures(inputs, n_workers=16)          # all CPUs by default
bad = [r for r in results if r["has_atomic_overlaps"]]

r = check_structure(path_or_structure_or_atoms)           # single structure

# subset to skip work: composition-only descriptors skip graph construction entirely
fast = check_structures(inputs, descriptors=["has_atomic_overlaps", "has_overcoordinated_c"])
```

Each result is a dict with `index`, `id`, `n_atoms`, and the requested descriptors. `DEFAULT_DESCRIPTORS` is the in-scope diagnostic suite (including bit-exact `has_high_charges`); `ALL_DESCRIPTORS` adds metadata, symmetry, and graph hashes. A structure that fails gets an `error` field (`on_error="record"`) instead of aborting the batch.

## ✅ Parity

Verified against a MOFChecker 2.0 checkout (used only as a behavioral oracle) via the harnesses in `scripts/`:

- **Real QMOFs:** 4500/4500 descriptor-comparisons (250 structures × 18 descriptors) — 100%.
- **Reference test CIFs:** 16/16.
- **Generated (distorted) structures:** 3899/3900. The single difference is **`has_lone_molecule`, where mofchecker-next is more correct** — see Limitations.

Reproduce: `scripts/qmof_parity.py` (real QMOFs), `scripts/generated_parity.py` (generated CIFs), `scripts/validate_subgraph_rx.py` (floating-solvent port).

## ⚠️ Limitations & deliberate differences

- **Healing not implemented.** `adding_hydrogen` / `adding_linker` raise `NotImplementedError`.
- **No porosity.** `is_porous` returns `None` (no bundled Zeo++).
- **`has_lone_molecule` is *more* correct than the reference.** MOFChecker 2.0's supercell+in-cell-filter heuristic silently misses finite molecules that wrap the unit-cell boundary (the origin-cell copy is truncated at the supercell face). mofchecker-next detects them via a topological finite-component test. This is the only descriptor that ever disagrees with the reference, only on pathological/distorted structures (0 disagreements on real QMOFs).
- **Graph construction is still the floor.** pymatgen's VESTA neighbor-finding (via `structuregraph_helpers`) is unchanged; the speedup is in the graph *algorithms*, not bond perception.
- **Determinism.** `symmetry_hash` is deterministic by design and will not match the reference's randomized value across runs.

## 🛠️ Install / build

```bash
python -m maturin develop --release          # build the Rust extension into the venv
python -m pytest -q                          # Python tests
cargo test --release --manifest-path rust/Cargo.toml   # Rust tests
```

Dependencies: `numpy`, `pymatgen` (loading + neighbor-finding), `structuregraph_helpers` (construction + hashes), `rustworkx`, `element-coder`, `libconeangle`.

## 🗂️ Layout

- `py/mofchecker_next/` — Python package (`checks/`, `diagnostics.py`, the `eqeq` subpackage).
- `py/mofchecker_next/checks/_subgraph_rx.py` — rustworkx floating-solvent + dimensionality.
- `rust/` — the `_rust` PyO3 extension (geometry + EQeq kernels).
- `scripts/` — parity harnesses and the speed benchmark.
- `tests/` — Rust and Python unit tests.
- `docs/DIAGNOSTIC_INVENTORY.md` — per-diagnostic parity status.

## ⚖️ Licensing

The **published package is GPLv2**, because it bundles `py/mofchecker_next/eqeq/` — a faithful translation of [EQeq](https://github.com/lsmo-epfl/EQeq) (GPLv2, see `py/mofchecker_next/eqeq/LICENSE`) — and the GPL governs the combined work. The non-eqeq sources are **MIT** (`LICENSE`); for an MIT-only build, omit the `eqeq` subpackage and the `has_high_charges` diagnostic.

The MOFChecker 2.0 checkout used as the behavioral oracle (ANCSA 1.0) is **not** redistributed; see `external/REFERENCE.md` to reproduce it locally.
