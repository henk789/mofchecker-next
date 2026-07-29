#!/usr/bin/env python3
"""Blinded adjudication of undercoordinated C/N disagreements between MOFChecker
0.9.6 and 2.0.

The two versions disagree on ``has_undercoordinated_c`` / ``has_undercoordinated_n``
for thousands of QMOF structures, and neither ships a labeled set, so "which is
more correct" cannot be settled from source alone. This harness turns the
disagreement into a hand-labelable experiment:

    sample   pick structures from each disagreement stratum (reproducible seed)
    build    extract the individual disagreeing ATOMS with their local geometry,
             shuffle them, and write a worksheet that does NOT say which version
             flagged which atom (provenance goes to key.json)
    score    read the filled verdicts and report, per version, how often its
             call matched the human label, plus a pool-weighted estimate and a
             sign test on the paired disagreements

Atom-level labeling is the point: "is this specific carbon missing a bond or
hydrogen?" is decidable by eye, while "is this structure valid?" is not.

Reference index files come from ``undercoordination_indices_ref.py``, run once
per pinned environment. Every stage is protocol-matched: both versions and this
script use the CIF's own site order (``symprec=None, primitive=False``), so atom
indices mean the same thing everywhere.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import random
from collections import Counter
from pathlib import Path

STRATA = (
    ("undercoordinated_c", "v096_only"),
    ("undercoordinated_c", "v20_only"),
    ("undercoordinated_n", "v096_only"),
    ("undercoordinated_n", "v20_only"),
)
VERDICTS = ("needs_bond", "ok", "unsure")


def _load(path: Path) -> dict[str, dict]:
    return {r["id"]: r for r in (json.loads(line) for line in path.open())}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


# -- sample ----------------------------------------------------------------
def sample(args) -> None:
    """Pick structures per stratum from paired boolean runs of both versions."""
    old, new = _load(args.old), _load(args.new)
    common = sorted(old.keys() & new.keys())
    pools: dict[str, list[str]] = {}
    for check, _ in STRATA[::2]:
        field = f"has_{check}"
        pools[f"{check}|v096_only"] = [
            i for i in common if old[i].get(field) is True and new[i].get(field) is False
        ]
        pools[f"{check}|v20_only"] = [
            i for i in common if new[i].get(field) is True and old[i].get(field) is False
        ]

    rng = random.Random(args.seed)
    selected: dict[str, list[str]] = {}
    for name, pool in pools.items():
        picks = sorted(pool) if len(pool) <= args.per_stratum else rng.sample(sorted(pool), args.per_stratum)
        selected[name] = sorted(picks)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "seed": args.seed,
        "per_stratum": args.per_stratum,
        "n_paired_structures": len(common),
        "pool_sizes": {name: len(pool) for name, pool in pools.items()},
        "selected": selected,
        "old_run": {"path": str(args.old), "sha256": _sha256(args.old)},
        "new_run": {"path": str(args.new), "sha256": _sha256(args.new)},
    }
    (args.out_dir / "sample.json").write_text(json.dumps(manifest, indent=2) + "\n")
    ids = sorted({i for picks in selected.values() for i in picks})
    (args.out_dir / "ids.txt").write_text("\n".join(ids) + "\n")
    print(json.dumps(manifest["pool_sizes"], indent=2))
    print(f"selected {len(ids)} structures -> {args.out_dir/'ids.txt'}")


# -- build -----------------------------------------------------------------
def _structure(path: Path):
    import warnings

    from pymatgen.core import Structure

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return Structure.from_file(str(path))


def _smallest_ring(simple_graph, index: int) -> int | None:
    import networkx as nx

    neighbors = list(simple_graph.neighbors(index))
    trimmed = simple_graph.copy()
    trimmed.remove_node(index)
    best = None
    for a, b in itertools.combinations(neighbors, 2):
        try:
            length = nx.shortest_path_length(trimmed, a, b)
        except (nx.NetworkXNoPath, nx.NodeNotFound):
            continue
        size = length + 2
        best = size if best is None else min(best, size)
    return best


def _atom_case(structure, graph, simple_graph, counts, index: int) -> dict:
    import numpy as np

    site = structure[index]
    neighbors = graph.get_connected_sites(index)
    entries = []
    for neighbor in neighbors:
        entries.append(
            {
                "element": str(neighbor.site.specie),
                "distance": round(float(structure.get_distance(index, neighbor.index)), 3),
                "neighbor_cn": counts[neighbor.index],
            }
        )
    angles, angles_from_distances = [], []
    for a, b in itertools.combinations(range(len(neighbors)), 2):
        v1 = site.coords - neighbors[a].site.coords
        v2 = site.coords - neighbors[b].site.coords
        cosine = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
        angles.append(round(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))), 1))
        # 2.0 derives this angle from three minimum-image distances instead of
        # the bonded images, which disagrees across cell boundaries. Show both so
        # a reviewer can tell a chemistry call from a PBC artifact.
        side_a = structure.get_distance(index, neighbors[a].index)
        side_b = structure.get_distance(index, neighbors[b].index)
        opposite = structure.get_distance(neighbors[a].index, neighbors[b].index)
        cosine_d = (side_a**2 + side_b**2 - opposite**2) / (2 * side_a * side_b)
        angles_from_distances.append(
            round(math.degrees(math.acos(max(-1.0, min(1.0, round(cosine_d, 6))))), 1)
        )
    return {
        "element": str(site.specie),
        "coordination_number": counts[index],
        "neighbors": entries,
        "angles_at_atom": sorted(angles),
        "angles_from_minimum_image_distances": sorted(angles_from_distances),
        "n_hydrogen_neighbors": sum(1 for e in entries if e["element"] == "H"),
        "n_metal_neighbors": sum(
            1 for neighbor in neighbors if bool(getattr(neighbor.site.specie, "is_metal", False))
        ),
        "smallest_ring": _smallest_ring(simple_graph, index),
    }


def _fragment_xyz(structure, index: int, radius: float) -> str:
    site = structure[index]
    lines = [f"{str(site.specie)} 0.000 0.000 0.000  # the atom under review"]
    for neighbor in structure.get_sites_in_sphere(site.coords, radius):
        offset = neighbor.coords - site.coords
        if float(abs(offset).sum()) < 1e-6:
            continue
        lines.append(f"{str(neighbor.specie)} {offset[0]:.3f} {offset[1]:.3f} {offset[2]:.3f}")
    return f"{len(lines)}\ncentered on reviewed atom\n" + "\n".join(lines) + "\n"


def build(args) -> None:
    import networkx as nx

    from mofchecker_next.checks.graph import build_structure_graph, connected_site_counts

    manifest = json.loads((args.review_dir / "sample.json").read_text())
    old, new = _load(args.old_indices), _load(args.new_indices)

    wanted: dict[tuple[str, str], set[str]] = {}
    for name, ids in manifest["selected"].items():
        check, stratum = name.split("|")
        wanted[(check, stratum)] = set(ids)

    cases = []
    skipped = Counter()
    for structure_id in sorted({i for ids in manifest["selected"].values() for i in ids}):
        old_record, new_record = old.get(structure_id), new.get(structure_id)
        if not (old_record and new_record and old_record.get("ok") and new_record.get("ok")):
            skipped["reference_error"] += 1
            continue
        structure = _structure(args.cif_dir / f"{structure_id}.cif")
        graph = build_structure_graph(structure)
        counts = connected_site_counts(graph)
        simple = nx.Graph()
        simple.add_nodes_from(range(len(structure)))
        simple.add_edges_from((int(u), int(v)) for u, v, _ in graph.graph.edges(data=True) if u != v)

        for check, stratum in STRATA:
            if structure_id not in wanted.get((check, stratum), ()):
                continue
            old_set, new_set = set(old_record[check]), set(new_record[check])
            only = sorted(old_set - new_set) if stratum == "v096_only" else sorted(new_set - old_set)
            if not only:
                skipped[f"no_atom_level_disagreement:{check}:{stratum}"] += 1
                continue
            for index in only[: args.atoms_per_structure]:
                case = {
                    "structure": structure_id,
                    "atom_index": index,
                    "check": check,
                    "flagged_by": "0.9.6" if stratum == "v096_only" else "2.0",
                    **_atom_case(structure, graph, simple, counts, index),
                }
                case["fragment"] = _fragment_xyz(structure, index, args.radius)
                cases.append(case)

    random.Random(manifest["seed"]).shuffle(cases)
    fragments = args.review_dir / "fragments"
    fragments.mkdir(parents=True, exist_ok=True)
    key, blinded = [], []
    for number, case in enumerate(cases, 1):
        case_id = f"case-{number:03d}"
        (fragments / f"{case_id}.xyz").write_text(case.pop("fragment"))
        key.append(
            {
                "case": case_id,
                "structure": case["structure"],
                "atom_index": case["atom_index"],
                "check": case["check"],
                "flagged_by": case.pop("flagged_by"),
            }
        )
        blinded.append({"case": case_id, **{k: v for k, v in case.items() if k != "structure"}})

    with (args.review_dir / "cases.jsonl").open("w") as handle:
        for case in blinded:
            handle.write(json.dumps(case) + "\n")
    (args.review_dir / "key.json").write_text(json.dumps(key, indent=2) + "\n")
    (args.review_dir / "worksheet.md").write_text(_worksheet(blinded))
    (args.review_dir / "verdicts.tsv").write_text(
        "case\tverdict\tnote\n" + "".join(f"{c['case']}\t\t\n" for c in blinded)
    )
    print(f"{len(blinded)} atom cases -> {args.review_dir}/worksheet.md (skipped: {dict(skipped)})")


def _worksheet(cases: list[dict]) -> str:
    lines = [
        "# Undercoordination adjudication worksheet",
        "",
        "Fill `verdict` in `verdicts.tsv` with one of "
        f"{', '.join(VERDICTS)} for every case, then run the `score` subcommand.",
        "",
        "`needs_bond` = this atom is genuinely missing a bond or hydrogen.",
        "`ok` = the local environment is chemically complete as drawn.",
        "`unsure` = excluded from scoring.",
        "",
        "Which checker version flagged each atom is deliberately not shown here;",
        "it lives in `key.json`. Geometry fragments are in `fragments/<case>.xyz`",
        "(reviewed atom first, at the origin).",
        "",
    ]
    for case in cases:
        neighbors = ", ".join(
            f"{n['element']} at {n['distance']} A (CN {n['neighbor_cn']})" for n in case["neighbors"]
        )
        lines += [
            f"## {case['case']}",
            "",
            f"- atom: {case['element']}, CN {case['coordination_number']}",
            f"- neighbors: {neighbors or 'none'}",
            f"- angles at atom (true bonded geometry): {case['angles_at_atom'] or 'n/a'}",
            "- same angles from minimum-image distances: "
            f"{case['angles_from_minimum_image_distances'] or 'n/a'}",
            f"- H neighbors: {case['n_hydrogen_neighbors']}, metal neighbors: {case['n_metal_neighbors']}",
            f"- smallest ring through atom: {case['smallest_ring'] or 'none'}",
            "",
        ]
    return "\n".join(lines)


# -- score -----------------------------------------------------------------
def _wilson(successes: int, total: int) -> tuple[float, float]:
    if total == 0:
        return (0.0, 1.0)
    z = 1.959963984540054
    phat = successes / total
    denominator = 1 + z**2 / total
    center = (phat + z**2 / (2 * total)) / denominator
    margin = z * math.sqrt(phat * (1 - phat) / total + z**2 / (4 * total**2)) / denominator
    return (max(0.0, center - margin), min(1.0, center + margin))


def _sign_test(favor_old: int, favor_new: int) -> float:
    """Two-sided exact binomial p-value for p=0.5 on decided pairs."""
    n = favor_old + favor_new
    if n == 0:
        return 1.0
    k = min(favor_old, favor_new)
    tail = sum(math.comb(n, i) for i in range(0, k + 1)) / 2**n
    return min(1.0, 2 * tail)


def score_review(key: list[dict], verdicts: dict[str, str], pool_sizes: dict[str, int], sampled: dict[str, int]) -> dict:
    """Score both versions against human labels on the adjudicated atoms."""
    per_version = {"0.9.6": Counter(), "2.0": Counter()}
    strata = Counter()
    favor = Counter()
    for entry in key:
        verdict = verdicts.get(entry["case"], "").strip()
        if verdict not in ("needs_bond", "ok"):
            strata["unlabeled_or_unsure"] += 1
            continue
        flagged_by = entry["flagged_by"]
        other = "2.0" if flagged_by == "0.9.6" else "0.9.6"
        truth_needs_bond = verdict == "needs_bond"
        # Only the flagging version calls this atom a problem, by construction.
        per_version[flagged_by]["tp" if truth_needs_bond else "fp"] += 1
        per_version[other]["fn" if truth_needs_bond else "tn"] += 1
        favor[flagged_by if truth_needs_bond else other] += 1
        stratum = f"{entry['check']}|{'v096_only' if flagged_by == '0.9.6' else 'v20_only'}"
        strata[f"{stratum}|{'real' if truth_needs_bond else 'spurious'}"] += 1

    result: dict = {"n_labeled": sum(favor.values()), "skipped": dict(strata), "versions": {}}
    for version, counts in per_version.items():
        correct = counts["tp"] + counts["tn"]
        total = sum(counts.values())
        low, high = _wilson(correct, total)
        result["versions"][version] = {
            **{k: counts[k] for k in ("tp", "fp", "fn", "tn")},
            "accuracy_on_disagreements": (correct / total) if total else 0.0,
            "accuracy_95ci": [round(low, 4), round(high, 4)],
        }
    result["favored_calls"] = {"0.9.6": favor["0.9.6"], "2.0": favor["2.0"]}
    result["sign_test_p"] = round(_sign_test(favor["0.9.6"], favor["2.0"]), 6)

    # Pool-weighted extrapolation: how many flags in each disagreement pool are real?
    weighted = {}
    for name, pool in pool_sizes.items():
        real = strata[f"{name}|real"]
        spurious = strata[f"{name}|spurious"]
        labeled = real + spurious
        if not labeled:
            continue
        rate = real / labeled
        low, high = _wilson(real, labeled)
        weighted[name] = {
            "pool": pool,
            "labeled": labeled,
            "real_rate": round(rate, 4),
            "real_rate_95ci": [round(low, 4), round(high, 4)],
            "estimated_real_flags_in_pool": round(pool * rate, 1),
            "sampled_structures": sampled.get(name, 0),
        }
    result["pool_weighted"] = weighted
    return result


def score(args) -> None:
    key = json.loads((args.review_dir / "key.json").read_text())
    manifest = json.loads((args.review_dir / "sample.json").read_text())
    verdicts = {}
    for line in (args.review_dir / "verdicts.tsv").read_text().splitlines()[1:]:
        if not line.strip():
            continue
        fields = line.split("\t")
        verdicts[fields[0].strip()] = fields[1] if len(fields) > 1 else ""
    result = score_review(
        key,
        verdicts,
        manifest["pool_sizes"],
        {name: len(ids) for name, ids in manifest["selected"].items()},
    )
    (args.review_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


def selfcheck(_args) -> None:
    key = [
        {"case": "case-001", "check": "undercoordinated_c", "flagged_by": "0.9.6"},
        {"case": "case-002", "check": "undercoordinated_c", "flagged_by": "0.9.6"},
        {"case": "case-003", "check": "undercoordinated_n", "flagged_by": "2.0"},
        {"case": "case-004", "check": "undercoordinated_n", "flagged_by": "0.9.6"},
    ]
    verdicts = {
        "case-001": "needs_bond",  # 0.9.6 right
        "case-002": "ok",  # 0.9.6 wrong
        "case-003": "needs_bond",  # 2.0 right
        "case-004": "unsure",  # excluded
    }
    pools = {"undercoordinated_c|v096_only": 1000, "undercoordinated_n|v20_only": 8}
    out = score_review(key, verdicts, pools, {})
    assert out["n_labeled"] == 3, out
    assert out["versions"]["0.9.6"] == {
        **{"tp": 1, "fp": 1, "fn": 1, "tn": 0},
        "accuracy_on_disagreements": 1 / 3,
        "accuracy_95ci": out["versions"]["0.9.6"]["accuracy_95ci"],
    }
    assert out["versions"]["2.0"]["tp"] == 1 and out["versions"]["2.0"]["fn"] == 1
    assert out["favored_calls"] == {"0.9.6": 1, "2.0": 2}
    assert out["skipped"]["unlabeled_or_unsure"] == 1
    weighted = out["pool_weighted"]["undercoordinated_c|v096_only"]
    assert weighted["labeled"] == 2 and weighted["real_rate"] == 0.5
    assert weighted["estimated_real_flags_in_pool"] == 500.0
    assert _sign_test(0, 0) == 1.0 and _sign_test(10, 0) < 0.01
    print("selfcheck ok")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("sample", help="pick structures from each disagreement stratum")
    p.add_argument("--old", type=Path, required=True, help="0.9.6 descriptor JSONL")
    p.add_argument("--new", type=Path, required=True, help="2.0 descriptor JSONL")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--per-stratum", type=int, default=15)
    p.add_argument("--seed", type=int, default=0)
    p.set_defaults(func=sample)

    p = sub.add_parser("build", help="extract blinded per-atom cases")
    p.add_argument("--review-dir", type=Path, required=True)
    p.add_argument("--old-indices", type=Path, required=True)
    p.add_argument("--new-indices", type=Path, required=True)
    p.add_argument("--cif-dir", type=Path, required=True)
    p.add_argument("--atoms-per-structure", type=int, default=1)
    p.add_argument("--radius", type=float, default=4.5)
    p.set_defaults(func=build)

    p = sub.add_parser("score", help="score filled verdicts")
    p.add_argument("--review-dir", type=Path, required=True)
    p.set_defaults(func=score)

    sub.add_parser("selfcheck", help="verify the scoring math").set_defaults(func=selfcheck)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
