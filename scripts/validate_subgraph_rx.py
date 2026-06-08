"""Validate the rustworkx finite-component port against the reference
structuregraph_helpers.get_subgraphs_as_molecules(return_unique=False) over many
structures. Parity = identical SET of molecule index-sets per structure.

  ~/projects/mofchecker-next/.venv/bin/python scripts/validate_subgraph_rx.py \
      --dirs DIR1 DIR2 --per-dir 150 --qmof-n 300
"""
from __future__ import annotations

import argparse
import json
import os
import random
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from qmof_parity import _repo_root  # noqa: E402

# Set QMOF_DIR to a local QMOF reference CIF directory (not redistributed here).
QMOF_DIR = Path(os.environ.get("QMOF_DIR", "qmof_cifs"))


def _as_set(list_of_lists):
    return frozenset(frozenset(int(i) for i in lst) for lst in list_of_lists)


def _one(path):
    from pymatgen.core import Structure

    from mofchecker_next.checks._subgraph_rx import finite_component_indices
    from mofchecker_next.checks.graph import build_structure_graph
    from structuregraph_helpers.subgraph import get_subgraphs_as_molecules

    s = Structure.from_file(path)
    sg = build_structure_graph(s, "vesta")
    ref_idx = get_subgraphs_as_molecules(sg, return_unique=False)[2]
    new_idx = finite_component_indices(sg)
    ref_set, new_set = _as_set(ref_idx), _as_set(new_idx)
    return {
        "id": Path(path).name,
        "match": ref_set == new_set,
        "bool_match": (len(ref_set) > 0) == (len(new_set) > 0),  # has_lone_molecule
        "ref_raw": len(ref_idx), "new_raw": len(new_idx),
        "ref_uniq": len(ref_set), "new_uniq": len(new_set),
        "only_ref": [sorted(x) for x in (ref_set - new_set)],
        "only_new": [sorted(x) for x in (new_set - ref_set)],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dirs", nargs="*", default=[])
    ap.add_argument("--per-dir", type=int, default=0)
    ap.add_argument("--qmof-n", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workers", type=int, default=12)
    args = ap.parse_args()

    paths = []
    for d in args.dirs:
        cifs = sorted(str(p) for p in Path(d).glob("*.cif"))
        paths.extend(cifs[: args.per_dir] if args.per_dir else cifs)
    if args.qmof_n:
        le150 = Path("/tmp/qmof_le150.json")
        pool = sorted(json.load(open(le150))) if le150.exists() else sorted(p.name for p in QMOF_DIR.glob("*.cif"))
        random.Random(args.seed).shuffle(pool)
        paths.extend(str(QMOF_DIR / nm) for nm in pool[: args.qmof_n])

    print(f"validating subgraph_rx on {len(paths)} structures, {args.workers} workers")
    from multiprocessing import Pool
    with Pool(args.workers) as pool:
        rows = pool.map(_one, paths)

    n = len(rows)
    matches = sum(r["match"] for r in rows)
    bool_matches = sum(r["bool_match"] for r in rows)
    set_mismatch = [r for r in rows if not r["match"]]
    bool_mismatch = [r for r in rows if not r["bool_match"]]
    print(f"\n  has_lone_molecule (DIAGNOSTIC) parity: {bool_matches}/{n} ({100*bool_matches/max(n,1):.3f}%)")
    print(f"  exact molecule-index-set parity:       {matches}/{n} ({100*matches/max(n,1):.3f}%)")
    for r in bool_mismatch[:25]:
        print(f"    BOOL FLIP {r['id']}: ref_uniq={r['ref_uniq']} new_uniq={r['new_uniq']}")
    if not bool_mismatch:
        print("  -> DIAGNOSTIC PARITY HOLDS (has_lone_molecule identical on all structures)")


if __name__ == "__main__":
    main()
