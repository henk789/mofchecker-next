#!/usr/bin/env python3
"""Dump per-atom undercoordinated C/N indices from an official MOFChecker env.

Both versions are constructed with ``symprec=None, angle_tolerance=None,
primitive=False`` so each keeps the CIF's own site order. Their ``from_cif``
defaults differ (0.9.6 primitive=False, 2.0 primitive=True), which would renumber
atoms and make cross-version index comparison meaningless.
"""

from __future__ import annotations

import argparse
import json
import warnings
from importlib.metadata import version
from pathlib import Path

warnings.filterwarnings("ignore")

CHECKS = {
    "undercoordinated_c": "no_undercoordinated_carbon",
    "undercoordinated_n": "no_undercoordinated_nitrogen",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cif-dir", type=Path, required=True)
    parser.add_argument("--ids", type=Path, required=True, help="one CIF stem per line")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    from mofchecker import MOFChecker
    from pymatgen.core import Structure

    def build(path: Path):
        kwargs = {
            "structure": Structure.from_file(path),
            "symprec": None,
            "angle_tolerance": None,
            "primitive": False,
        }
        try:
            return MOFChecker(**kwargs)
        except TypeError as exc:
            if "linker_structure" not in str(exc):
                raise
            return MOFChecker(**kwargs, linker_structure=None)

    versions = {
        "mofchecker": version("mofchecker"),
        "pymatgen": version("pymatgen"),
        "structuregraph_helpers": version("structuregraph-helpers"),
    }
    protocol = {"symprec": None, "angle_tolerance": None, "primitive": False}
    stems = [line.strip() for line in args.ids.read_text().split() if line.strip()]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as output:
        for stem in stems:
            record: dict = {
                "id": stem,
                "ok": False,
                "versions": versions,
                "protocol": protocol,
            }
            try:
                checker = build(args.cif_dir / f"{stem}.cif")
                record["n_sites"] = len(checker.structure)
                for name, key in CHECKS.items():
                    record[name] = sorted(int(i) for i in checker.checks[key].flagged_indices)
                record["ok"] = True
            except Exception as exc:  # noqa: BLE001
                record["error"] = f"{type(exc).__name__}: {exc}"[:300]
            output.write(json.dumps(record) + "\n")
            output.flush()
    print(f"wrote {len(stems)} records to {args.out}")


if __name__ == "__main__":
    main()
