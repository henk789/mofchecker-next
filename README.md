# mofchecker-next

Incremental Rust/Python replacement kernels for the MOFChecker MOF diagnostics.

Python owns CIF/structure loading, pymatgen integration, diagnostics
orchestration, and parity scripts. Rust owns small deterministic kernels
(minimum-image distances, short contacts, neighbor candidates, connected
components, graph degrees, and the EQeq charge-equilibration solve).

## Status

All in-scope checker diagnostics reach parity with MOFChecker 2.0, verified by
the parity scripts under `scripts/`:

- Atomic overlaps, over/under-coordination (C/N/H, rare-earth, alkali/alkaline),
  floating solvent, 3D-connectivity, suspicious terminal oxo, geometrically
  exposed metal, simple composition descriptors.
- `has_high_charges` via a faithful Rust port of EQeq (bit-exact charges).
- `possible_charged_fused_ring`, `positive_charge_from_linkers`,
  `negative_charge_from_linkers`, and open metal sites (`has_oms`).

Validation: 16/16 on the MOFChecker reference CIFs and 9000/9000
descriptor-comparisons (500 QMOF structures x 18 descriptors) at 100%.

## Drop-in `MOFChecker`

A MOFChecker-compatible class exposes the same properties and
`get_mof_descriptors` API, so existing code can switch over:

```python
from mofchecker_next import MOFChecker

mc = MOFChecker.from_cif("structure.cif")     # also .from_ase(atoms) / MOFChecker(structure)
mc.has_atomic_overlaps, mc.has_oms, mc.metal_number
mc.graph_hash, mc.spacegroup_symbol
descriptors = mc.get_mof_descriptors()        # OrderedDict of all descriptors
```

Implemented: all diagnostics (with `*_indices` accessors), composition, metadata
(`formula`, `density`), symmetry (`spacegroup_symbol/number`, `symmetry_hash`),
and graph hashes (`graph_hash`, `undecorated_graph_hash`, `scaffold_hash`,
`undecorated_scaffold_hash`). The graph is built once and reused across checks.

Differences from the reference: `adding_hydrogen`/`adding_linker` (healing) raise
`NotImplementedError`; `is_porous` returns `None` (no bundled Zeo++);
`symmetry_hash` is deterministic (the reference's depends on Python hash
randomization and is not reproducible across runs).

## Batch validation

`mofchecker_next.batch` runs the full diagnostic set over many structures
efficiently — useful for validating generated structures (e.g. from a diffusion
model). It accepts pymatgen `Structure`, ASE `Atoms`, or CIF paths (mixed),
builds each structure's graph once, and parallelizes across structures with
`multiprocessing`.

```python
from mofchecker_next.batch import check_structures, check_structure

# inputs may be paths, pymatgen Structures, ASE Atoms, or a mix
results = check_structures(inputs, n_workers=16)          # all CPUs by default
overlapping = [r for r in results if r["has_atomic_overlaps"]]

# one structure
r = check_structure(atoms_or_structure_or_path)

# pick a subset (dropping has_oms removes the most expensive check;
# a composition-only subset skips graph construction entirely)
fast = check_structures(inputs, descriptors=["has_atomic_overlaps", "has_overcoordinated_c"])
```

Each result is a dict with `index`, `id`, `n_atoms`, and the requested
descriptors. `DEFAULT_DESCRIPTORS` is the full in-scope suite, including the
bit-exact EQeq `has_high_charges`; `ALL_DESCRIPTORS` additionally returns the
explicit `oms_indices` list. Failed structures get an `error` field
(`on_error="record"`) instead of aborting the batch.

## Layout

- `py/mofchecker_next/` — Python package (checks, diagnostics, the `eqeq` subpackage).
- `rust/` — the `_rust` PyO3 extension (geometry + EQeq kernels).
- `scripts/` — local parity harnesses against the reference oracle.
- `tests/` — Rust and Python unit tests.
- `docs/DIAGNOSTIC_INVENTORY.md` — per-diagnostic parity status.

## Build & test

```bash
python -m maturin develop --release   # build the Rust extension into the venv
python -m pytest -q                   # Python tests
cargo test --release --manifest-path rust/Cargo.toml   # Rust tests
```

## Licensing

This project is **MIT-licensed** (see `LICENSE`), with one exception: the
subpackage `py/mofchecker_next/eqeq/` is a faithful translation of
[EQeq](https://github.com/lsmo-epfl/EQeq) and vendors its data tables. EQeq is
**GPLv2**, so that subpackage (and any distribution bundling it) is GPLv2 — see
`py/mofchecker_next/eqeq/LICENSE`. For an MIT-only distribution, omit the `eqeq`
subpackage and the `has_high_charges` diagnostic.

The MOFChecker 2.0 checkout used as a behavioral oracle (ANCSA 1.0) is **not**
redistributed here; see `external/REFERENCE.md` to reproduce it locally.
