from __future__ import annotations

import argparse
import json
from pathlib import Path

from mofchecker_next.batch import DEFAULT_DESCRIPTORS, check_structures, summarize_results


def main() -> None:
    p = argparse.ArgumentParser(description="Check MOF CIFs with persistent workers.")
    p.add_argument("paths", nargs="*", type=Path, help="CIF files or directories containing *.cif")
    p.add_argument("--cif_dir", type=Path, help="Directory containing *.cif")
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--n_workers", type=int, default=1)
    p.add_argument("--chunksize", type=int, default=1)
    p.add_argument("--timeout_s", type=float, default=120.0, help="per-CIF timeout; 0 disables")
    p.add_argument("--output_json", type=Path)
    p.add_argument("--all_descriptors", action="store_true")
    p.add_argument("--no_progress", action="store_true")
    args = p.parse_args()

    inputs: list[Path] = []
    roots = ([args.cif_dir] if args.cif_dir else []) + args.paths
    for root in roots:
        if root.is_dir():
            inputs.extend(sorted(root.glob("*.cif")))
        else:
            inputs.append(root)
    if args.limit is not None:
        inputs = inputs[: args.limit]
    if not inputs:
        raise SystemExit("no CIFs found")

    descriptors = None if args.all_descriptors else list(DEFAULT_DESCRIPTORS)
    results = check_structures(
        inputs,
        n_workers=args.n_workers,
        descriptors=descriptors,
        chunksize=args.chunksize,
        progress=not args.no_progress,
        timeout_s=args.timeout_s or None,
    )
    summary = summarize_results(results)
    print(json.dumps(summary, indent=2, sort_keys=True))
    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
