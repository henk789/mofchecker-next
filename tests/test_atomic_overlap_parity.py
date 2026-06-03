import json
import subprocess
import sys
from pathlib import Path

import pytest
from pymatgen.core import Lattice, Structure

from mofchecker_next.checks.geometry import build_overlap_cutoff_matrix, check_atomic_overlaps


def involved_indices(diagnostics):
    return sorted({atom.index for diagnostic in diagnostics for atom in diagnostic.atoms})


def test_artificial_positive_overlap_matches_reference_boolean_and_indices(tmp_path):
    root = Path(__file__).resolve().parents[1]
    if not (root / ".venv-ref" / "bin" / "python").exists():
        pytest.skip("MOFChecker 2.0 reference environment .venv-ref is not available")

    structure = Structure(
        Lattice.cubic(10.0),
        ["C", "C"],
        [[0.0, 0.0, 0.0], [0.01, 0.0, 0.0]],
    )
    cif_path = tmp_path / "positive_overlap.cif"
    structure.to(filename=str(cif_path))

    out_path = tmp_path / "reference.json"
    subprocess.run(
        [sys.executable, str(root / "scripts" / "run_reference_mofchecker.py"), str(cif_path), "--out", str(out_path)],
        check=True,
    )
    reference = json.loads(out_path.read_text())["reference"]

    cutoff_matrix = build_overlap_cutoff_matrix([6, 6], {6: 0.76})
    diagnostics = check_atomic_overlaps(structure, cutoff_matrix)

    assert reference["has_atomic_overlaps"] is True
    assert len(diagnostics) > 0
    assert sorted(reference["overlapping_indices"]) == involved_indices(diagnostics)
    assert diagnostics[0].atoms[0].image == (0, 0, 0)
    assert diagnostics[0].atoms[1].image == (0, 0, 0)
    assert diagnostics[0].values["distance_angstrom"] < diagnostics[0].values["cutoff_angstrom"]
