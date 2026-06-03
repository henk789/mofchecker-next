# Diagnostic Inventory

This inventory tracks checker diagnostics targeted for parity with MOFChecker 2.0. It excludes healing/correction outputs.

## Implemented / Scaffolded

| Diagnostic | Reference API | Reference Source | Reference Output | Current Status | Parity Target |
|---|---|---|---|---|---|
| Atomic overlaps | `has_atomic_overlaps`, `get_overlapping_indices()` | `checks/local_structure/overlapping_atoms.py` | boolean + flagged atom indices | Rust-backed short-contact diagnostics implemented | boolean + involved atom index set |
| Lone molecule / floating solvent | `has_lone_molecule`, `lone_molecule_indices` | `checks/floating_solvent.py` | boolean + index lists | production wrapper uses `structuregraph_helpers.subgraph`; local parity reporting added | boolean + flattened/index-list parity |
| Overcoordinated C | `has_overcoordinated_c`, `overvalent_c_indices` | `checks/local_structure/overcoordinated_carbon.py` | boolean + flagged C indices | production wrapper uses `structuregraph_helpers` graph plus Rust degree helper | flagged index set |
| Overcoordinated N | `has_overcoordinated_n`, check flagged indices | `checks/local_structure/overcoordinated_nitrogen.py` | boolean + flagged N indices | production wrapper uses `structuregraph_helpers` graph plus Rust degree helper | flagged index set |
| Overcoordinated H | `has_overcoordinated_h`, `overvalent_h_indices` | `checks/local_structure/overcoordinated_hydrogen.py` | boolean + flagged H indices | pymatgen neighbor-count helper implemented | flagged index set |

## Simple Global Diagnostics

| Diagnostic | Reference API | Reference Source | Reference Output | Current Status | Parity Target |
|---|---|---|---|---|---|
| Has carbon | `has_carbon` | `checks/global_structure/__init__.py`, `checks/utils/get_indices.py` | boolean | Python helper implemented; local parity reporting added | boolean |
| Has hydrogen | `has_hydrogen` | `checks/global_structure/__init__.py`, `checks/utils/get_indices.py` | boolean | Python helper implemented; local parity reporting added | boolean |
| Has nitrogen | `has_nitrogen` | `checks/global_structure/__init__.py`, `checks/utils/get_indices.py` | boolean | Python helper implemented; local parity reporting added | boolean |
| Has metal | `has_metal` | `checks/global_structure/__init__.py`, `definitions.py` | boolean | Python helper implemented with caller-supplied metal set; local parity reporting added | boolean |
| Metal number | `metal_number` | `mofchecker/__init__.py`, `checks/utils/get_indices.py` | integer | Python helper implemented with MOFChecker's metals-minus-Sb convention; local parity reporting added | integer |

## Charge Diagnostics

| Diagnostic | Reference API | Reference Source | Reference Output | Current Status | Parity Target |
|---|---|---|---|---|---|
| High charges | `has_high_charges` | `checks/charge_check.py` (EQeq via `pyeqeq`) | boolean (`None` if EQeq missing) | implemented as the GPL `mofchecker_next.eqeq` subpackage: a faithful Rust translation of EQeq's C++ + vendored ionization/charge-center tables | boolean + bit-exact per-atom charges |

Notes:

- The `eqeq` subpackage is **GPLv2** (EQeq is GPL); the rest of the project is MIT. See `py/mofchecker_next/eqeq/LICENSE`.
- `scripts/compare_high_charges_parity.py` verifies bit-exact charge parity (max|Δq| = 0.000000) against the `pyeqeq` oracle across all local reference CIFs.
- EQeq's hand-rolled CIF parser reads single-letter elements from the label column as hydrogen; MOFChecker 2.0's `has_high_charges` inherits this bug because it round-trips through a pymatgen CIF. Our kernel reads the correct elements, so the parity script feeds the oracle an EQeq-safe CIF (label = element symbol). The diagnostic boolean matches the reference 16/16 regardless.

## Graph / Dimensionality Diagnostics

| Diagnostic | Reference API | Reference Source | Reference Output | Current Status | Parity Target |
|---|---|---|---|---|---|
| 3D connected graph | `has_3d_connected_graph` | `checks/global_structure/graphcheck.py` | boolean | Python wrapper uses pymatgen dimensionality on `structuregraph_helpers` graph; local parity reporting added | boolean |

## Higher-Risk Local Structure Diagnostics

These require branch-by-branch inspection before implementation.

| Diagnostic | Reference API | Reference Source | Current Status |
|---|---|---|---|
| Undercoordinated C | `has_undercoordinated_c`, `undercoordinated_c_indices` | `checks/local_structure/undercoordinated_carbon.py` | index parity implementation added; candidate H positions remain out of diagnostic scope |
| Undercoordinated N | `has_undercoordinated_n`, `undercoordinated_n_indices` | `checks/local_structure/undercoordinated_nitrogen.py` | index parity implementation added; candidate H positions remain out of diagnostic scope |
| Undercoordinated rare earth | `has_undercoordinated_rare_earth`, `undercoordinated_rare_earth_indices` | `checks/local_structure/undercoordinated_rare_earth.py` | parity helper added |
| Undercoordinated alkali/alkaline | `has_undercoordinated_alkali_alkaline` | `checks/local_structure/undercoordinated_alkaline.py` | parity helper added |
| Geometrically exposed metal | `has_geometrically_exposed_metal`, `geometrically_exposed_metal_indice` | `checks/local_structure/geometrically_exposed_metal.py` | parity helper and local CIF reporting added using `element-coder` + `libconeangle` |
| Suspicious terminal oxo | `has_suspicious_terminal_oxo`, `suspicious_terminal_oxo_indices` | `checks/local_structure/false_oxo.py` | parity helper and local CIF reporting added; reference-visible indices are metals, not terminal O atoms |
| Possible charged fused ring | `possible_charged_fused_ring` | `checks/local_structure/fused_ring.py` | implemented in `checks/charge_oms.py`; boolean parity verified (16/16 reference CIFs + synthetic benzimidazole positive case) |
| Positive linker charge | `positive_charge_from_linkers` | `checks/local_structure/positive_charge.py` | implemented in `checks/charge_oms.py`; int parity verified (16/16 + synthetic Ge/Sb, quaternary-N, O-CN3 triggers) |
| Negative linker charge | `negative_charge_from_linkers` | `checks/local_structure/negative_charge.py` | implemented in `checks/charge_oms.py`; int parity verified (15/16 reference CIFs carry non-zero values exercising O/S/P/halogen/N branches) |
| OMS | `has_oms`, `oms_indice` | `checks/oms/` | implemented in `checks/charge_oms.py`; boolean + index parity verified (5 open / 1 closed-metal / 10 no-metal on reference CIFs). No-metal cells return no-OMS instead of raising NoMetal |

Parity for these four is checked by `scripts/compare_charge_oms_parity.py` against the reference oracle. They are faithful ports of the reference graph heuristics, running on the same `structuregraph_helpers` "vesta" graph.

## Explicitly Out Of Scope For Now

| Feature | Reason |
|---|---|
| `adding_hydrogen` | healing/correction output with CIF-writing side effects |
| `adding_linker` | healing/correction output with CIF-writing side effects |
| corrected CIF writing | outside diagnostic-only scope |
| porosity / Zeo++ | optional external binary dependency |
| graph/scaffold/symmetry hashes | descriptors, not checker diagnostics |
| GUI/report generation | outside Rust kernel and diagnostics parity scope |
