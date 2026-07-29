# Plan for an actually corrected `next` mode

## Decision

Keep `0.9.6` and `2.0` as immutable compatibility profiles. Add
`mode="next"` as an explicitly experimental corrected profile; every report must
also record a resolved schema such as `next-dev-1` plus the source SHA. Once its
validation gates pass, freeze it as `next-1.0`. Published results use the frozen
name, never the rolling `next` alias.

Do not make `next` the default while it is being designed. Mofgen remains on
exact 0.9.6 for literature/model-selection continuity and evaluates `next` in
shadow on the same saved CIFs until `next-1.0` is frozen.

## Historical shadow baseline: `next-dev-1`

The metrics below describe the superseded broad compatibility-like composite and
remain archived for provenance. Rolling `mode="next"` now resolves to
`next-dev-2`, whose conservative named composite excludes unvalidated
undercoordination and other application-dependent report-only checks. Production
gates and implementation status are maintained in `NEXT_PRODUCTION_PLAN.md`.

Implemented in the historical baseline:

- immutable mode profiles with `next` resolving to provisional `next-dev-1`;
- common non-mutating input path for CIF/Structure/ASE;
- corrected periodic finite components and parallel-edge coordination counts;
- bonded-image CN2-carbon angles (2.0's minimum-image-distance behavior remains
  isolated to explicit compatibility mode);
- true-element, cached Rust EQeq with explicit `total_charge=0`, JSON-safe
  per-atom charges, max/sum/residual/threshold metrics, support errors, and
  fail-closed charge-sum validation;
- per-structure mode/profile metadata and valid/invalid/indeterminate status;
- batch summaries with unconditional validity, error categories, and
  per-descriptor error counts;
- mofgen `--mofchecker_mode {0.9.6,2.0,next}` while retaining 0.9.6 as default.

Initial enriched-256 QMOF smoke audit: 0.9.6 remains at 71/256 = 27.734% with
0/4,096 legacy parity disagreements; explicit 2.0 is 206/256 = 80.469%;
`next-dev-1` is 208/256 = 81.250%. Both modern modes had zero checker errors; the
next difference is two corrected CN2-carbon PBC calls. Rust EQeq scored all 256
with zero errors, maximum charge-sum residual 4.27e-14, and no `|q| > 4` cases.

Generated-CIF shadow audit on identical saved samples preserves the established
ranking and stays numerically close to the modern metric:

| saved 1,000-CIF set | 0.9.6 unconditional | 2.0 unconditional | next-dev-1 unconditional | next indeterminate |
|---|---:|---:|---:|---:|
| medium inter (`3802456`) | 49.3% | 56.6% | 56.6% | 14 |
| large inter (`3802457`) | 52.6% | 57.6% | 57.6% | 13 |

Every next indeterminate case is an explicit unsupported-actinide Rust EQeq
result (Pa/U/Np/Pu); none silently passed. Conditional next rates are 57.40% and
58.36%, respectively. Thus large still beats medium under both unconditional and
conditional scoring. A separate older medium-none CIF set also gave identical
2.0/next unconditional validity (33.6%, eight unsupported-charge indeterminates).

Full-QMOF `next-dev-1` audit over all 20,374 saved CIFs: 17,954 valid, 2,234
invalid, and 186 indeterminate, giving 88.122% unconditional and 88.934%
conditional validity. Every composite error is an unsupported true-element EQeq
case: U 169, Th 10, Pu 4, Np 3. The supported 20,188 structures had no
`max |q| > 4` failures; median/p95/p99/max were 1.207/2.398/3.215/3.359 and the
maximum charge-sum residual was 2.64e-13. At threshold 3, 304 supported ground
truth structures would fail, so retaining 4 is currently the less destructive
provisional choice. Generated medium/large inter samples still show 1.22%/0.71%
high-charge flags at threshold 4, so the corrected charge check is active rather
than vacuously passing.

These distributions remain useful monitoring evidence, not release calibration.
The remaining production gates are listed in `NEXT_PRODUCTION_PLAN.md`.

## First principle: stop pretending one Boolean answers five questions

The current composite mixes chemistry, dataset membership, and application
preferences. `next` reports separate statuses:

1. `structure_status`: parseable cell, finite coordinates, supported occupancy,
   no impossible overlap, internally consistent periodic graph.
2. `local_chemistry_status`: plausible local valence/coordination for checks that
   are applicable and validated.
3. `charge_status`: balanced, inconsistent, ambiguous, unsupported, or not run.
4. `component_status`: framework, stray atom, finite molecule/solvent, and
   dimensionality.
5. `scope_status`: whether it belongs to a declared application domain, e.g.
   carbon-containing, metal-containing, desolvated generated MOFs.

Porosity, 3D connectivity, open metal sites, and hydrogen presence are reported
properties, not universal validity conditions. Legitimate MOFs can be 1D/2D,
nonporous, H-free, solvated, or contain intended open metal sites.

A scalar is allowed only as a named application profile, initially
`generated_desolvated_organic_mof_v1`, whose exact required fields are stored in
the output. Generic `valid=True` is removed from corrected mode.

## Phase 0 — freeze history and establish the evaluation contract

**Work**

- Preserve the existing full-QMOF 0.9.6 parity fixture: 20,374 CIFs, 16 fields,
  325,984/325,984 matches, 16,259 valid.
- Preserve explicit 2.0 behavior separately; do not let corrected fixes leak into
  it unless they are required for 2.0 parity.
- Add `next` to the mode registry, but initially return only the new report
  schema—no corrected scalar yet.
- Replace scattered `if mode == ...` policy decisions with one small immutable
  profile record containing input normalization, element sets, check functions,
  descriptor set, and composite name.
- Every batch summary records mode, resolved profile, package version/SHA,
  input protocol, composite fields, requested/scored/error counts, and error
  categories.

**Gate**

- 0.9.6 remains at zero full-QMOF disagreements.
- Explicit 2.0 regression fixtures remain unchanged.
- Every output can be reproduced from its manifest without inferring defaults.

## Phase 1 — trustworthy inputs and periodic geometry

**Corrected `next` policy**

- CIF, pymatgen, and ASE inputs all enter the same normalization path.
- Preserve the submitted cell and atom indices: no implicit primitive conversion,
  symmetrization, atom deletion, or coordinate snapping.
- Reject/mark indeterminate malformed lattices, NaN/Inf coordinates, unsupported
  disorder/partial occupancy, and nonsensical site species.
- Keep periodic-image edge identity and parallel edges. Coordination number is
  the number of bonded images, not unique atom indices.
- Compute angles from the actual bonded image vectors. Never reconstruct an
  angle from three independently minimum-imaged distances.
- Make every edge expose element pair, distance, image offset, and bond-rule
  provenance. Mark near-cutoff bonds ambiguous instead of presenting graph
  inference as exact chemistry.

**Tests**

- Invariance under atom permutation, integer cell translations, origin shifts,
  equivalent supercells, and CIF/Structure/ASE round trips.
- Focused periodic-image, parallel-edge, wrapped-fragment, and near-cutoff cases.
- Synthetic malformed/disordered inputs must return explicit indeterminate
  results, not fabricated Booleans or crashes.

**Gate**: 100% invariance and expected-result agreement on these mechanical tests.

## Phase 2 — high-confidence structural checks

Implement corrected checks as structured results with
`status = pass|fail|indeterminate|not_applicable`, severity, atom indices,
periodic neighbors, measured values, rule/threshold, and explanation.

Start with checks whose defect has an objective geometric definition:

- atomic overlaps;
- impossible H overcoordination;
- egregious C/N overcoordination;
- finite disconnected components, split into stray atoms and multi-atom
  molecules;
- framework dimensionality as a property, not a failure;
- cell/coordinate/occupancy sanity.

Do **not** initially put undercoordinated C/N, the *legacy parser-based* charge
Boolean, exposed metal, terminal oxo, alkali/rare-earth coordination, OMS,
porosity, or hydrogen presence into the corrected scalar. Correct true-element
Rust EQeq lands in Phase 4 and is required before `next-1.0`; the others are
chemistry- or application-dependent and need separate evidence.

**Validation corpus**

- curated known-good QMOF and BWDB structures;
- saved generated CIFs from mofgen;
- controlled corruptions with known targets: duplicate/move/delete one atom,
  stretch one bond, inject a finite fragment, remove a counterion, distort the
  cell, and create a wrapped component;
- preserve pristine and corrupted pairs so each check is tested against the
  exact defect it claims to detect.

**Gate**

- Mechanical checks: lower 95% confidence bound for precision and recall ≥0.95
  on held-out controlled defects.
- No tuning on the held-out test split.
- Report per-check confusion matrices; aggregate validity alone is insufficient.

## Phase 3 — rebuild C/N undercoordination from labels, not inherited code

The existing 53-case blinded disagreement worksheet is a pilot, not enough to
freeze a rule.

**Expand it**

- 0.9.6-only and 2.0-only C/N atoms;
- agreement-positive and agreement-negative controls;
- pristine structures plus atom-deletion corruptions;
- at least 200–300 atom cases, stratified by CN, metal neighbor, ring membership,
  nitrile/alkyne motif, and periodic-boundary involvement;
- two independent chemistry reviewers, blinded to checker/version/corruption,
  followed by adjudication; record inter-rater agreement.

**Implement only after labeling**

- distinguish covalent neighbors from metal coordination;
- recognize nitrile/alkyne and aromatic/pyridine-like environments;
- handle charged/deprotonated sites explicitly rather than guessing from angle;
- keep Boolean status, flagged indices, explanation, and any repair position
  generated from one result object, so they cannot disagree;
- treat ambiguous graph/bond-order cases as indeterminate.

**Gate**

- inter-rater Cohen's κ ≥0.8 before treating labels as ground truth;
- held-out lower 95% confidence bound: precision ≥0.90 and recall ≥0.80;
- if either gate fails, undercoordination remains reported but stays outside the
  corrected scalar. No threshold is chosen because it raises aggregate validity.

## Phase 4 — fix and use the Rust EQeq charge check

The solver itself is not the problem. The in-tree Rust EQeq implementation is a
faithful translation of EQeq and already matches `pyeqeq` on all 16/16 reference
CIFs when both receive the same correctly identified elements (worst
`max |Δq| = 2.22e-16`, far below EQeq's 3-decimal output precision). The published MOFChecker
bug is the CIF-label input path around that solver.

### Corrected `next` implementation

- Call the Rust kernel directly from the pymatgen site's true element symbol;
  never write/reparse a CIF and never enable `reference_cif_labels`.
- Keep the verified EQeq parameters, Ewald settings, and Rust numerics unchanged.
- Require an explicit expected total charge in the charge API. The
  `generated_desolvated_organic_mof_v1` profile declares `total_charge=0`; a
  generic profile with unknown total charge returns indeterminate rather than
  silently assuming neutrality.
- Return every per-atom partial charge, max absolute charge, sum of charges,
  expected total charge, solver settings, and high-charge atom indices.
- Add/return solver residual or at minimum verify `sum(q)` against the requested
  total charge within a pinned tolerance.
- Elements absent from the vendored EQeq table return `unsupported`, never the
  Z=0/hydrogen fallback. Extend the table only from a cited, versioned parameter
  source with dedicated numerical tests.
- Name the result `eqeq_charge_plausibility`, not formal charge or charge balance:
  EQeq returns partial charges under a supplied total-charge constraint.

### Calibrate the criterion instead of inheriting 3 or 4

The corrected check must be used in `next`, but its cutoff must earn its value.
Store `max_abs_eqeq_charge` continuously, then select and freeze the
`next-1.0` threshold using:

- curated pristine QMOF/BWDB structures;
- generated malformed structures;
- paired controlled corruptions, including atom/linker/counterion removal,
  severe bond stretching, and overlap creation;
- an element-stratified train/calibration/test split so a metal family cannot
  appear only on one side.

Choose the threshold on the calibration split before opening the held-out test.
Report the complete precision/recall curve and element/support coverage. Start
with one global threshold; add element-specific thresholds only if prespecified
and clearly supported, not to rescue individual failures.

**Gate for inclusion in `next-1.0`**

- true-element Rust-vs-`pyeqeq` numerical parity passes on all supported test
  structures;
- no silent fallback and zero unexplained solver failures;
- lower 95% confidence bound for precision ≥0.90 on held-out charge-plausibility
  defects, with recall and unsupported coverage reported;
- total-charge assumption is explicit in every result;
- the frozen threshold and EQeq settings are part of the resolved profile.

This corrected EQeq check is required by the mofgen generated-structure profile.
An indeterminate/unsupported charge result therefore makes that profile
indeterminate and fail closed in its unconditional validity denominator.

### Separate future formal charge-balance diagnostic

EQeq charge plausibility still does not prove formal oxidation/linker charge
balance because its total charge is an input. A later, separately named
`formal_charge_consistency` check may combine oxidation-state feasibility,
bond-valence evidence, and linker formal-charge inference. It must be validated
against curated oxidation-state/counterion labels and must return
`balanced | inconsistent | ambiguous | unsupported`. It supplements corrected
EQeq; it does not replace or mislabel it.

## Phase 5 — define and freeze `next-1.0`

Pre-register the first application composite before running final benchmarks.
A conservative starting profile for mofgen is:

### `generated_desolvated_organic_mof_v1`

Required:

- structure sanity passes;
- declared carbon-containing and metal-containing scope passes;
- no atomic overlap;
- no impossible H/C/N overcoordination;
- no stray atom or finite molecule (because this profile explicitly assumes
  desolvated generated structures);
- corrected true-element Rust EQeq charge plausibility passes under the profile's
  explicit `total_charge=0`, using the frozen calibrated threshold;
- validated C/N undercoordination passes, **only if Phase 3 passes its gate**.

Reported but not universally required:

- hydrogen presence;
- 1D/2D/3D connectivity;
- porosity;
- open/exposed metal sites;
- terminal oxo motifs;
- rare-earth/alkali coordination;
- formal charge-balance status until a separately validated method exists.

If a required check is indeterminate, the scalar is indeterminate and counts as
invalid in the unconditional denominator. It never silently passes.

Freeze the profile field list, thresholds, element sets, graph method, and input
normalization under `next-1.0`; subsequent chemistry changes become `next-1.1`
or `next-2.0`, never an in-place metric change.

## Phase 6 — mofgen rollout

1. Keep 0.9.6 as the selection/default metric.
2. Add an explicit mofgen `--mofchecker_mode next` and run it in shadow on the
   **same saved CIFs** as 0.9.6—never regenerate to compare evaluators.
3. Log separate namespaces: `mofchecker_v096/*` and
   `mofchecker_next_dev_1/*`, including resolved profile and error counts.
4. Compare model rankings across checkpoints, seeds, QMOF ground truth, BWDB
   ground truth, and generated CIFs.
5. Switch model selection only after `next-1.0` is frozen, its benchmark report is
   reviewed, and ranking stability is understood. Historical runs remain labeled
   and untouched.

## Likely code changes

- `py/mofchecker_next/core.py`: route through a profile instead of scattered mode
  branches; preserve Boolean compatibility properties.
- `py/mofchecker_next/profiles.py` (new, small): immutable definitions for
  `0.9.6`, `2.0`, rolling development profiles, and named composites.
- `py/mofchecker_next/diagnostics.py`: add structured check status/evidence and a
  report container while retaining existing `Diagnostic` records.
- `py/mofchecker_next/checks/graph.py`: corrected graph evidence, image-vector
  angles, ambiguity, and valence results.
- `py/mofchecker_next/eqeq/`: true-element Rust EQeq result with explicit total
  charge, support/error state, solver evidence, and calibrated plausibility flag.
- `py/mofchecker_next/batch.py` and `cli.py`: mode/profile metadata, tri-state
  summaries, unconditional denominator, and per-check errors.
- `tests/`: immutable parity tests, invariance/property tests, controlled
  corruption tests, and frozen labeled benchmark fixtures.

## Explicit non-goals

- Do not maximize validity percentage.
- Do not treat DFT-relaxed QMOF as automatically valid ground truth.
- Do not use Crystalite validity as MOFChecker ground truth.
- Do not silently repair structures before scoring them.
- Do not infer correctness from parity with either flawed historical version.
- Do not replace the verified Rust EQeq kernel with another dependency. Do not
  add a formal-charge model without a held-out benchmark and frozen artifact/version.
- Do not collapse ambiguous chemistry into `False` merely to preserve a Boolean
  API; compatibility modes already serve that API.

## Recommended first implementation slice

Implement only Phase 0 and Phase 1 plus the new structured report. That creates a
safe `mode="next"` skeleton with corrected input/periodic semantics and no false
claim of chemical validity. Then build the controlled-corruption benchmark and
expand the blinded C/N labels before adding chemistry to the composite. The
corrected Rust EQeq charge path and calibration follow before freezing
`next-1.0`; a separate formal charge-balance method can come later.
