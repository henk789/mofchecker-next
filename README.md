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
