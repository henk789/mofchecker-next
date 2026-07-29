# What a corrected MOFChecker should contain

Compatibility and correctness are separate products. `mode="0.9.6"` reproduces
the literature evaluator exactly, including its dependency/parser defects;
`mode="2.0"` targets the newer implementation. A future corrected mode must get
a new name/version and must not silently change either historical score.

## Required design

1. **Three explicit evaluator identities**
   - exact 0.9.6 compatibility;
   - exact 2.0 compatibility;
   - a separately versioned corrected composite.
   Every result records package version/SHA, mode, composite fields, input type,
   primitive/symmetry settings, denominator, and checker errors.

2. **One immutable composite per identity**
   Presence and problem fields are named and versioned. Adding or removing a
   field creates a new composite instead of changing old scores in place.

3. **Tri-state descriptors**
   Each check returns pass, fail, or indeterminate. Errors and unsupported
   chemistry are indeterminate and fail closed in the unconditional validity
   rate; missing values never pass silently.

4. **Correct periodic graph semantics**
   Coordination number counts bonded periodic-image edges, including parallel
   edges to different images of the same indexed site. Angles use the actual
   bonded image vectors, not three independently minimum-imaged distances.

5. **Validated C/N valence checks**
   The rules distinguish covalent from metal-coordination neighbors, recognize
   nitriles/alkynes and charged/deprotonated sites, and keep every Boolean flag,
   flagged index, and proposed repair position consistent. No dead CN branches
   or disabled `len(indices) == len(positions)` invariant. Thresholds are fitted
   or audited on a blinded atom-level set, not selected from aggregate validity.

6. **A corrected charge check using the verified Rust EQeq kernel**
   The in-tree solver is bit-exact with `pyeqeq` when given identical correct
   inputs, so keep it and fix the surrounding contract: pass true element
   identities directly, require/record expected total charge, expose per-atom
   charges and solver evidence, and return unsupported rather than a silent Z=0
   fallback. Calibrate and freeze the high-partial-charge threshold on held-out
   pristine/corrupted MOFs instead of inheriting 3 or 4. Call this EQeq charge
   plausibility—not formal charge balance, which requires a separate diagnostic.

7. **Explicit element sets**
   Metal, alkali/alkaline, and rare-earth sets are versioned data. A corrected
   rare-earth check should state whether it follows IUPAC (including Sc/Y) rather
   than inheriting a changing pymatgen convenience property.

8. **Correct finite-component detection**
   Stray atoms and multi-atom molecules are separate diagnostics computed with
   periodic image-offset topology, so components crossing a cell boundary are
   not lost. A legacy union can remain only in compatibility mode.

9. **Pure validity path**
   Evaluating the composite does not invoke healing, linker decomposition,
   porosity, hashes, or unrelated descriptors. Each descriptor is isolated, and
   one optional failure cannot invalidate unrelated checks.

10. **Trust-boundary robustness**
    CIF/ASE/pymatgen inputs follow a documented, identical normalization path;
    occupancy/disorder, malformed cells, unsupported elements, timeouts, memory
    exhaustion, and solver failures produce explicit errors rather than crashes
    or fabricated Booleans.

11. **Atom-level explanations**
    Every failure exposes atom indices, periodic neighbors, distances/angles,
    the rule and threshold that fired, and any uncertainty. This makes a score
    auditable instead of an opaque Boolean.

12. **Evidence before release**
    Tests include exact historical parity sets, periodic-boundary and
    parallel-edge regressions, unsupported-element charge cases, generated
    distorted CIFs, and a blinded expert-labeled C/N disagreement benchmark.
    Report per-descriptor precision/recall and uncertainty, not only a higher
    aggregate validity percentage.

## Proposed corrected composite

Keep the mechanically reliable presence, overlap, overcoordination,
finite-component, terminal-oxo, and exposed-metal checks after their periodic
and element-set fixes. Include C/N undercoordination only after the blinded audit
establishes useful precision/recall. Include corrected true-element Rust EQeq
charge plausibility only with an explicit total-charge assumption and a
held-out-calibrated threshold; unsupported/errors are indeterminate and fail
closed. Connectivity and porosity should remain reported properties, not
validity requirements, unless a separately named application profile explicitly
requires them.
