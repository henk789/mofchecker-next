#!/usr/bin/env python3
"""Reproducible paired benchmark for the objective atomic-overlap check.

Pristine CIFs are negative controls. Each paired corruption places one site at
another site's coordinate and is therefore an unambiguous positive. Splits are
assigned by pristine CIF bytes before corruption, preventing pair leakage.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from mofchecker_next.batch import check_structure
from mofchecker_next.core import normalize_structure
from mofchecker_next.profiles import profile_manifest


def wilson_lower(successes: int, total: int, z: float = 1.959963984540054) -> float:
    if total == 0:
        return 0.0
    p = successes / total
    denominator = 1 + z * z / total
    center = p + z * z / (2 * total)
    spread = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total))
    return (center - spread) / denominator


def split_for(data: bytes, seed: int) -> str:
    value = int.from_bytes(hashlib.sha256(str(seed).encode() + b"\0" + data).digest()[:8], "big") % 100
    return "train" if value < 70 else "calibration" if value < 85 else "test"


def overlap_corruption(structure):
    corrupted = structure.copy()
    if len(corrupted) == 1:
        corrupted.append(corrupted[0].species, corrupted[0].frac_coords)
    else:
        corrupted.replace(len(corrupted) - 1, corrupted[-1].species, corrupted[0].frac_coords)
    return corrupted


def confusion(labels: list[bool], predictions: list[bool | None]) -> dict:
    tp = sum(label and prediction is True for label, prediction in zip(labels, predictions))
    fp = sum(not label and prediction is True for label, prediction in zip(labels, predictions))
    tn = sum(not label and prediction is False for label, prediction in zip(labels, predictions))
    fn = sum(label and prediction is False for label, prediction in zip(labels, predictions))
    indeterminate = sum(prediction is None for prediction in predictions)
    precision_denominator = tp + fp
    recall_denominator = tp + fn
    return {
        "tp": tp, "fp": fp, "tn": tn, "fn": fn, "indeterminate": indeterminate,
        "precision": tp / precision_denominator if precision_denominator else 0.0,
        "recall": tp / recall_denominator if recall_denominator else 0.0,
        "precision_wilson95_lower": wilson_lower(tp, precision_denominator),
        "recall_wilson95_lower": wilson_lower(tp, recall_denominator),
    }


def collect_paths(roots: list[Path]) -> list[Path]:
    paths = []
    for root in roots:
        if root.is_dir():
            paths.extend(sorted(root.glob("*.cif")))
        elif root.suffix.lower() == ".cif":
            paths.append(root)
        else:
            paths.extend(Path(line.strip()) for line in root.read_text().splitlines() if line.strip() and not line.startswith("#"))
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--seed", type=int, default=20260729)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    paths = collect_paths(args.paths)
    if args.limit is not None:
        paths = paths[: args.limit]
    records = []
    grouped = {split: ([], []) for split in ("train", "calibration", "test")}
    for path in paths:
        data = path.read_bytes()
        split = split_for(data, args.seed)
        structure = normalize_structure(path, corrected=True)
        for corruption, label, candidate in (
            ("pristine", False, structure),
            ("duplicate_coordinate", True, overlap_corruption(structure)),
        ):
            result = check_structure(candidate, mode="next", descriptors=["has_atomic_overlaps"])
            prediction = result.get("has_atomic_overlaps")
            if not isinstance(prediction, bool):
                prediction = None
            grouped[split][0].append(label)
            grouped[split][1].append(prediction)
            records.append({
                "id": path.name,
                "cif_sha256": hashlib.sha256(data).hexdigest(),
                "split": split,
                "corruption": corruption,
                "expected_has_atomic_overlaps": label,
                "predicted_has_atomic_overlaps": prediction,
                "error": result.get("error") or result.get("errors", {}).get("has_atomic_overlaps"),
            })

    report = {
        "benchmark": "paired_atomic_overlap_v1",
        "seed": args.seed,
        "n_pristine": len(paths),
        "profile_manifest": profile_manifest("next"),
        "metrics": {split: confusion(*values) for split, values in grouped.items()},
        "records": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report["metrics"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
