# MOFChecker 2.0 Diagnostics Parity Plan

This plan targets checker diagnostics only. It does not target GUI, correction/healing workflows, CIF parsing in Rust, charge balancing, adding hydrogens/linkers, or report generation.

The MOFChecker 2.0 checkout under `external/mofchecker_2_ref/` is a behavioral oracle. Do not edit files there. Because its license is ANCSA 1.0, do not copy source blocks or fixtures into this project unless license compatibility is explicitly approved.

## Parity Definition

For a supported checker diagnostic, parity means:

- The reference and `mofchecker-next` return the same boolean result.
- Flagged atom index sets match when MOFChecker 2.0 exposes indices.
- Exceptions are classified as reference exception, our exception, unsupported, or known mismatch.
- Our diagnostics may include richer information, such as pair distances or periodic image vectors, but they must reduce to the reference-visible result.

## Architecture Rules

- Python owns CIF loading, pymatgen/ASE integration, graph construction, chemistry orchestration, diagnostics schema, CLI/reporting, and parity scripts.
- Rust owns small deterministic kernels only: MIC, periodic distances, short contacts, neighbor candidates, connected components, and simple graph degree operations.
- Rust must receive arrays, edges, atomic numbers, lattice matrices, and numeric cutoffs. It must not parse CIFs or infer chemistry.
- Production dependencies require approval before adding. Local parity scripts may use `.venv-ref` and the reference checkout.

## Current Completed Foundation

- Rust/Python MIC with `r_ij` convention: `delta_frac = frac_j + image_j - frac_i`.
- Rust/Python short-contact detection with caller-supplied cutoff matrix.
- Python atomic-overlap diagnostics.
- Atomic-overlap local parity script.
- Python cutoff matrix builder from caller-supplied radii.
- Rust/Python neighbor-candidate kernel.
- Rust/Python connected components over explicit edges.
- Rust/Python node degrees over explicit edges.
- Explicit-edge non-periodic component diagnostic scaffolding.
- Generic overcoordination-by-degree helper with caller-supplied exclusions.
- H overcoordination helper using pymatgen neighbor counting and caller-supplied H VDW radius.
- Suspicious terminal oxo helper using graph connectivity and the reference disallowed-metal list.
- Geometrically exposed metal helper using the reference cone-angle dependencies (`element-coder`, `libconeangle`).
- Local graph parity script using reference-built graph arrays.
- EQeq `has_high_charges` diagnostic: the `mofchecker_next.eqeq` subpackage is a faithful, bit-exact Rust translation of EQeq's C++ (`GetJ`/`Qeq`/`RoundCharges`) plus vendored ionization/charge-center tables. Verified max|Δq| = 0.000000 vs the `pyeqeq` oracle across all local reference CIFs. This subpackage is GPLv2 (see its LICENSE); the rest of the project remains MIT.

## Checkpoint 1: Diagnostic Inventory

Goal: define exactly which MOFChecker 2.0 diagnostics are in scope and how parity will be measured.

Tasks:

- Create a checker diagnostic inventory table.
- For each checker, record reference source files, public API property, output type, index accessor if available, dependencies, and parity target.
- Mark unsupported/deferred checks clearly.

Acceptance:

- Every diagnostics-only descriptor in scope has an inventory row.
- Healing/correction outputs are explicitly out of scope.

## Checkpoint 2: Robust Parity Harness

Goal: make local parity runs repeatable and machine-readable.

Tasks:

- Extend parity scripts to accept descriptor groups.
- Add `--fail-on-mismatch` where useful.
- Emit normalized JSON rows with reference values, our values, index sets, match status, and notes.
- Keep reference CIFs local under `external/`; do not copy them into committed fixtures.

Acceptance:

- One command can run local parity over `external/mofchecker_2_ref/test_cases/**/*.cif` for supported diagnostics.
- Mismatches are actionable.

## Checkpoint 3: Atomic Overlap Full Diagnostic Parity

Goal: complete the first production-quality diagnostic.

Tasks:

- Keep cutoff construction in Python.
- Compare against MOFChecker 2.0 `has_atomic_overlaps` and `get_overlapping_indices()`.
- Add artificial committed tests for same-cell, boundary, skew-cell, H, metal/nonmetal, and no-overlap cases.
- Keep image vectors as richer diagnostics with our `image_j` convention.

Acceptance:

- Boolean and index-set parity across local reference CIFs that run.
- Artificial positive/negative tests pass.

## Checkpoint 4: Graph Construction Decision

Goal: decide how production Python obtains graph edges for graph-based diagnostics.

Options:

- Use `structuregraph_helpers` in production Python for MOFChecker 2.0 parity.
- Keep `structuregraph_helpers` reference-only and implement a separate graph builder later.
- Use pymatgen graph builders directly and accept documented non-parity until proven.

Decision:

- Use `structuregraph-helpers` as a production Python dependency for diagnostics parity first. Package metadata reports MIT license for `structuregraph-helpers==0.0.9`.
- Evaluate whether selected parts should later be ported to Rust for speed after parity tests are stable.

Acceptance:

- Decision recorded.
- Production dependency added.
- Future speed-port evaluation remains explicit and separate from parity implementation.

### Future Speed-Port Evaluation For `structuregraph_helpers`

Candidate areas to evaluate after parity:

- edge extraction and graph degree operations, already partly covered by Rust kernels;
- connected components over explicit edges, already covered by Rust kernels;
- subgraph/floating-solvent supercell expansion and boundary filtering;
- local-environment graph construction via the `vesta` method.

Do not port graph construction until we have edge-set, `to_jimage`, connected-site, and descriptor parity tests proving equivalent behavior.

## Checkpoint 5: Overcoordinated H/C/N Parity

Goal: match MOFChecker 2.0 overcoordination diagnostics.

Reference behavior:

- C: graph coordination number `> 4`, excluding atoms with any metal neighbor or any boron neighbor.
- N: graph coordination number `> 4`, excluding atoms with any metal neighbor.
- H: pymatgen `Structure.get_neighbors(site, vdw_radius("H"))`; flag if neighbor count `> 1`.

Tasks:

- Build production wrappers after graph construction decision.
- Compare booleans and flagged index sets.
- Add artificial tests for C/N positive cases, metal exclusions, boron exclusion for C, H positive/negative cases.

Acceptance:

- Boolean and flagged-index parity across local reference CIFs.

## Checkpoint 6: Floating Solvent / Lone Molecule Parity

Goal: match `has_lone_molecule` and `lone_molecule_indices`.

Reference behavior:

- Uses `structuregraph_helpers.subgraph.get_subgraphs_as_molecules`.
- Expands to a 3x3x3 supercell.
- Finds connected components in an undirected graph.
- Rejects components with nonzero `to_jimage` boundary-crossing edges.
- Returns index lists for non-periodic components.

Tasks:

- Initially call `structuregraph_helpers` from Python if approved as a production dependency.
- Otherwise keep local parity via reference-built graph arrays and mark production implementation incomplete.
- Add artificial tests for isolated atom, isolated molecule, periodic component, and boundary-crossing component.

Acceptance:

- Boolean and flagged-index parity across local reference CIFs.

## Checkpoint 7: Simple Global Diagnostics

Goal: implement/check simple diagnostics that do not need Rust.

Targets:

- `has_carbon`
- `has_hydrogen`
- `has_nitrogen`
- `has_metal`
- `metal_number`

Tasks:

- Inspect reference index helpers and metal definitions.
- Implement Python helpers.
- Compare descriptor parity.

Acceptance:

- Parity across local reference CIFs.

## Checkpoint 8: 3D Connectedness

Goal: decide whether to expose `has_3d_connected_graph` parity.

Reference behavior:

- Uses pymatgen `get_dimensionality_larsen(structure_graph) == 3`.

Tasks:

- Keep Python-only initially.
- Requires graph construction decision.

Acceptance:

- Boolean parity across local reference CIFs if included in diagnostics scope.

## Checkpoint 9: Higher-Risk Local Structure Diagnostics

Order:

1. Undercoordinated carbon
2. Undercoordinated nitrogen
3. Undercoordinated rare earth
4. Undercoordinated alkali/alkaline
5. False terminal oxo
6. Geometrically exposed metal
7. Fused ring
8. Positive/negative linker charge diagnostics
9. OMS

For each diagnostic:

- Inspect source.
- Document constants, tolerances, inputs, outputs, and side effects.
- Add branch-level artificial tests.
- Add local CIF parity comparison.
- Implement minimal Python/Rust pieces only where deterministic kernels help.

Notes:

- False terminal oxo is implemented. MOFChecker 2.0 reports the metal index for this check, not the terminal oxygen index.
- Geometrically exposed metal is implemented after dependency approval. The current `libconeangle` package requires `float64` arrays, so inputs are normalized before calling the same cone-angle API.

Acceptance:

- One diagnostic reaches parity before starting the next.

Status: Checkpoint 9 complete. All nine higher-risk diagnostics reach parity. Fused ring, positive/negative linker charge, and OMS are implemented in `checks/charge_oms.py` as faithful ports of the reference graph heuristics (running on the same `structuregraph_helpers` "vesta" graph) and verified by `scripts/compare_charge_oms_parity.py` (16/16 reference CIFs plus synthetic positive-direction triggers). With `has_high_charges` also done, every in-scope checker diagnostic now matches the reference.

## Deferred Out Of Scope

- `adding_hydrogen`
- `adding_linker`
- writing corrected CIFs
- charge balancing/correction workflows
- porosity / Zeo++ unless explicitly requested
- GUI
- broad CLI/report generation
- graph/hash/symmetry descriptors unless needed for diagnostics parity

## Current Next Step

Start Checkpoint 1 by creating a concrete diagnostic inventory for all diagnostics currently considered in scope.
