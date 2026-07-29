# Adjudicating undercoordinated C/N: 0.9.6 vs 2.0

MOFChecker 0.9.6 and 2.0 disagree about `has_undercoordinated_c` and
`has_undercoordinated_n` on thousands of QMOF structures, and 2.0's higher
validity rate is driven mostly by those two checks. Neither release ships a
labeled set, and the 2.0 paper's manual audit is not a version ablation, so
source reading alone cannot say which behavior is right. This harness decides it
empirically with a small blinded hand-labeling study.

Tooling: `scripts/undercoordination_review.py` (`sample`, `build`, `score`,
`selfcheck`) plus `scripts/undercoordination_indices_ref.py`, which runs inside a
pinned reference environment. Scoring math is covered by
`tests/test_undercoordination_review.py`.

## Design

- **Unit of judgment is the atom, not the structure.** "Is this carbon missing a
  hydrogen?" is decidable by eye; "is this MOF valid?" is not.
- **Protocol-matched.** Both references are constructed with `symprec=None`,
  `angle_tolerance=None`, `primitive=False`, so both keep the CIF's own site
  order and indices are comparable. Their `from_cif` defaults differ (0.9.6
  `primitive=False`, 2.0 `primitive=True`) and would otherwise renumber atoms.
- **Blinded.** The worksheet never says which version flagged an atom; provenance
  lives in `key.json`. This matters because the expected answer is guessable from
  a version's known tendencies.
- **Stratified.** Four strata: {C, N} x {flagged only by 0.9.6, flagged only by
  2.0}. The 2.0-only strata are tiny, so they are sampled to exhaustion and the
  0.9.6-only strata are subsampled with a recorded seed.
- **Both geometry readings are shown.** The angle from the bonded periodic images
  (physically correct) and the angle 2.0 reconstructs from three minimum-image
  distances. When they differ, the disagreement is a PBC artifact rather than a
  chemistry judgment, and the reviewer can see that.

## Rubric

For each case, decide whether the reviewed atom is missing a bond or hydrogen.

`needs_bond` — the environment is incomplete, for example:

- aromatic or sp2 carbon with two ring neighbors near 1.39 A and an angle near
  120 degrees, with no hydrogen: a missing ring H;
- sp3 carbon or nitrogen with fewer neighbors than its geometry implies;
- amine nitrogen with one H where the geometry implies two.

`ok` — the environment is complete as drawn, for example:

- nitrile or alkyne carbon, roughly linear, with a short (about 1.15-1.25 A)
  triple bond;
- carboxylate, amide, or similar carbon with three neighbors;
- pyridine-type ring nitrogen with two ring bonds and no hydrogen;
- deprotonated N or O whose charge is a charge-check matter, not a missing bond;
- a bent M-C or M-N coordination bond, which is soft and need not be linear.

`unsure` — anything you cannot call; excluded from scoring rather than guessed.

## Running it

```bash
REVIEW=$AUDIT/review/undercoordination

# 1. sample structures from each disagreement stratum
python scripts/undercoordination_review.py sample \
  --old  $AUDIT/results/qmof_full_20374_composite/v096_zatom.jsonl \
  --new  $AUDIT/results/qmof_full_20374_composite/v20_zatom.jsonl \
  --out-dir $REVIEW --per-stratum 15 --seed 0

# 2. per-atom flagged indices from each pinned reference environment
$AUDIT/envs/v096/bin/python scripts/undercoordination_indices_ref.py \
  --cif-dir $QMOF_DIR --ids $REVIEW/ids.txt --out $REVIEW/indices_v096.jsonl
$AUDIT/envs/v20/bin/python  scripts/undercoordination_indices_ref.py \
  --cif-dir $QMOF_DIR --ids $REVIEW/ids.txt --out $REVIEW/indices_v20.jsonl

# 3. blinded worksheet + xyz fragments
python scripts/undercoordination_review.py build --review-dir $REVIEW \
  --old-indices $REVIEW/indices_v096.jsonl \
  --new-indices $REVIEW/indices_v20.jsonl --cif-dir $QMOF_DIR

# 4. fill verdicts.tsv by hand, then
python scripts/undercoordination_review.py score --review-dir $REVIEW
```

## What the score means

`score` reports, per version, true/false positives and negatives over the
adjudicated atoms, accuracy with a Wilson interval, and an exact sign test on the
paired disagreements. Because only disagreement atoms are sampled, these are
accuracies *conditional on the versions disagreeing*, not overall accuracies.

`pool_weighted` extrapolates each stratum's labeled "real problem" rate back to
its full pool, estimating how many of that version's exclusive flags across all
of QMOF are genuine. That is the number that says whether 2.0's validity gain is
better chemistry or under-detection.

## Current instance

Sampled from all 20,374 released QMOF CIFs under the matched protocol. Pool
sizes: 1,040 C and 1,315 N flagged only by 0.9.6; 41 C and 8 N flagged only by
2.0. That asymmetry is itself the headline: 2.0 almost never sees a problem
0.9.6 misses, so the study is mostly asking how many of 0.9.6's extra flags are
real. 53 atom cases are pending labels.
