# `next` production plan

## Release decision

Do not rename the rolling profile to `next-1.0` until every required scientific
gate below has evidence. Engineering work proceeds as `next-dev-2`; exact
`0.9.6` and explicit `2.0` remain separate compatibility identities and the
default remains `0.9.6`.

The production composite is named
`generated_desolvated_organic_mof_v1`. Its required fields are frozen in each
result instead of inferred from the rolling alias.

## Plan and status

| Gate | Requirement | Status |
|---|---|---|
| Compatibility lock | Immutable 0.9.6/2.0 profile fixtures; no corrected behavior selected by the rolling alias name | Implemented; full 0.9.6 oracle remains 0/325,984 disagreements, post-refactor enriched fixture remains 0/4,096 |
| Honest composite | Require input sanity, C/metal scope, overlap, egregious C/N/H overcoordination, finite components, and corrected EQeq only | Implemented in `next-dev-2`; inherited C/N undercoordination, H presence, exposed metal, terminal oxo, rare-earth/alkali rules, OMS and connectivity are report-only |
| Trust boundary | Reject malformed/non-finite/singular/disordered/partial-occupancy/non-element inputs without fabricated Booleans | Implemented and tested; corrected CIF parsing disables coordinate snapping and primitive conversion |
| Periodic invariance | Atom permutation, integer translations, origin shifts, equivalent supercells and CIF round trips preserve mechanical labels | Implemented regression tests; integer-translation component bug found and fixed by a wrapped graph-only computational view that preserves submitted cell/index order |
| Structured evidence | Named check status, severity, rule, atoms/images, values, error/unsupported state; named composite derived from required checks | Implemented for every corrected-composite requirement; `valid` remains only a documented compatibility alias |
| Reproducibility | Package version, Git SHA/dirty state, source-tree digest, exact profile, settings, input protocol, descriptors, composite fields and denominators | Implemented in per-result metadata and batch manifest; CLI can emit ordered result JSONL and hashes ordered input paths plus bytes |
| Error accounting | Missing/non-Boolean/errored required fields are indeterminate; report-only errors are inert; duplicate basenames cannot overwrite errors | Implemented and tested |
| Objective mechanical benchmark | Frozen paired corruptions, split before corruption, per-check confusion matrices and Wilson 95% bounds | Atomic-overlap v1 implemented. On 1,000 QMOF CIFs its held-out test split is 154/154 TP and 154/154 TN, precision/recall 1.0, lower 95% bounds 0.9757, zero indeterminate. Additional component/overcoordination corpora remain required before claiming their empirical gates. |
| EQeq implementation | True elements, explicit total charge, no fallback, residual validation, per-atom evidence, numerical parity | Implemented; 16/16 `pyeqeq` parity, worst max-charge deviation 2.22e-16 |
| EQeq scientific calibration | Element-stratified pristine/corrupted train/calibration/test corpus; threshold selected before held-out test; precision lower bound >=0.90 | **Blocked on labeled charge-plausibility corpus.** Threshold 4 remains provisional monitoring policy, not calibrated evidence. Unsupported actinides remain explicit `unsupported`/indeterminate. |
| C/N undercoordination | 200–300 blinded atom labels, two reviewers, kappa >=0.8, held-out precision/recall gates | **Blocked on human chemistry labels.** The inherited checks remain report-only and cannot affect the production composite. |
| Release freeze | Clean source tree, immutable `next-1.0`, wheel/sdist smoke, all manifests/hashes, same-CIF QMOF/BWDB/generated shadow review | Wheel/sdist build and clean-wheel install smoke passed locally when source-build libclang prerequisites were supplied; release CI now repeats installed-wheel manifest/Rust smoke tests. Freeze remains blocked by mechanical and EQeq gates. C/N labels are optional because C/N undercoordination is excluded. |

## Shadow metric after the composite correction

The same saved CIFs were rescored; these values are not substituted for missing
per-check validation:

| set | valid | invalid | indeterminate | unconditional | conditional |
|---|---:|---:|---:|---:|---:|
| full QMOF (20,374) | 18,164 | 2,024 | 186 | 89.153% | 89.974% |
| medium inter generated (1,000) | 660 | 326 | 14 | 66.0% | 66.94% |
| large inter generated (1,000) | 666 | 321 | 13 | 66.6% | 67.48% |

Every composite indeterminate remains an explicit unsupported-actinide EQeq
case. The large model still ranks above the medium model. The increase from the
historical `next-dev-1` composite is expected: unvalidated undercoordination,
exposed-metal and application-dependent chemistry rules no longer silently
control the corrected scalar.

## Required work before `next-1.0`

1. Curate and freeze controlled finite-component and C/N/H-overcoordination pairs;
   require held-out lower 95% precision and recall bounds at least 0.95 for each
   mechanical check retained in the composite.
2. Build the charge-plausibility corpus and freeze one global EQeq threshold only
   after calibration and held-out evaluation. Keep unsupported elements
   indeterminate unless a cited parameter source is added and numerically tested.
3. Re-run exact full-QMOF 0.9.6 parity from the final clean release candidate.
4. Add immutable `mode="next-1.0"`; point rolling `next` at the same profile only
   after benchmark review. Keep `0.9.6` as default.
5. Build wheel and sdist from a clean tree, install each in an empty environment,
   run profile/invariance/charge tests, and archive the source-tree and benchmark
   hashes.
6. Score identical saved full-QMOF, BWDB and generated CIFs under 0.9.6, 2.0 and
   next-1.0; report unconditional/conditional denominators, every unsupported or
   error category, and model-ranking changes.

## Explicit release blockers

No code can manufacture the missing charge labels or independent chemistry
reviews. Until the charge gate passes, `next-dev-2` is suitable for shadow
production evaluation and engineering validation, but publishing it as a frozen
scientific `next-1.0` metric would be misleading.
