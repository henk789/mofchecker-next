from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from mofchecker_next.batch import (
    DEFAULT_DESCRIPTORS,
    LEGACY_DEFAULT_DESCRIPTORS,
    NEXT_DEFAULT_DESCRIPTORS,
    _stream_worker_from_env,
    check_cif_paths,
    summarize_results,
)
from mofchecker_next.core import DEFAULT_MODE, MODES


def main() -> None:
    p = argparse.ArgumentParser(description="Check MOF CIFs with persistent workers.")
    p.add_argument("paths", nargs="*", type=Path, help="CIF files, list files, or directories containing *.cif")
    p.add_argument("--cif_dir", type=Path, help="Directory containing *.cif")
    p.add_argument("--input_list", type=Path, help="Text file with one CIF path per line")
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--n_workers", type=int, default=1)
    p.add_argument("--chunksize", type=int, default=1)
    p.add_argument("--timeout_s", type=float, default=120.0, help="hard per-CIF timeout; 0 disables")
    p.add_argument("--output_json", type=Path, help="aggregate summary JSON")
    p.add_argument("--output_results_jsonl", type=Path, help="optional ordered per-structure results")
    p.add_argument("--all_descriptors", action="store_true")
    p.add_argument("--mode", choices=MODES, default=DEFAULT_MODE, help="evaluation profile (default: 0.9.6)")
    p.add_argument("--no_progress", action="store_true")
    p.add_argument("--soft_timeout", action="store_true", help="legacy in-process alarm instead of kill/respawn workers")
    p.add_argument("--_stream_worker", action="store_true", help=argparse.SUPPRESS)
    args = p.parse_args()
    if args._stream_worker:
        _stream_worker_from_env()
        return

    inputs: list[Path] = []
    roots = ([args.cif_dir] if args.cif_dir else []) + args.paths
    if args.input_list:
        roots += [Path(line.strip()) for line in args.input_list.read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")]
    for root in roots:
        if root.is_dir():
            inputs.extend(sorted(root.glob("*.cif")))
        elif root.suffix == ".cif":
            inputs.append(root)
        else:
            inputs.extend(Path(line.strip()) for line in root.read_text().splitlines() if line.strip() and not line.lstrip().startswith("#"))
    if args.limit is not None:
        inputs = inputs[: args.limit]
    if not inputs:
        raise SystemExit("no CIFs found")

    defaults = {
        "0.9.6": LEGACY_DEFAULT_DESCRIPTORS,
        "2.0": DEFAULT_DESCRIPTORS,
        "next": NEXT_DEFAULT_DESCRIPTORS,
    }[args.mode]
    descriptors = None if args.all_descriptors else list(defaults)
    results = check_cif_paths(
        inputs,
        n_workers=args.n_workers,
        descriptors=descriptors,
        progress=not args.no_progress,
        timeout_s=args.timeout_s or None,
        hard_timeout=not args.soft_timeout,
        mode=args.mode,
    )
    summary = summarize_results(results)
    summary["mode"] = args.mode
    digest = hashlib.sha256()
    for path in inputs:
        digest.update(str(path.resolve()).encode() + b"\0")
        digest.update(hashlib.sha256(path.read_bytes()).digest())
    summary["input_manifest"] = {
        "n_inputs": len(inputs),
        "ordered_paths_and_bytes_sha256": digest.hexdigest(),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    if args.output_results_jsonl:
        args.output_results_jsonl.parent.mkdir(parents=True, exist_ok=True)
        with args.output_results_jsonl.open("w") as handle:
            for result in results:
                handle.write(json.dumps(result, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
