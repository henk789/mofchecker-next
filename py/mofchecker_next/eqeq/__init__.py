"""Clean-room EQeq-style charge equilibration subpackage.

Python owns the chemistry (element parameters, charge centers, orchestration);
the Rust ``_rust.eqeq_charges`` kernel owns the deterministic numerics (Ewald
summation and the linear solve). Used by the ``has_high_charges`` diagnostic to
match MOFChecker 2.0's behaviour without depending on the GPL-licensed EQeq/pyeqeq.
"""

from __future__ import annotations

import numpy as np

from mofchecker_next.diagnostics import AtomRef, Diagnostic

from .parameters import (
    DEFAULT_H_ELECTRON_AFFINITY,
    DEFAULT_LAMBDA,
    AtomParameters,
    parameters_for,
)

__all__ = [
    "compute_charges",
    "has_high_charges",
    "high_charge_indices",
    "check_high_charges",
    "DEFAULT_LAMBDA",
    "DEFAULT_H_ELECTRON_AFFINITY",
    "DEFAULT_CHARGE_THRESHOLD",
]

# MOFChecker 2.0 flags a structure when any |charge| exceeds this threshold.
DEFAULT_CHARGE_THRESHOLD = 4.0

# EQeq method defaults (main.cpp / pyeqeq CLI).
DEFAULT_ETA = 50.0
DEFAULT_REAL_CELLS = 2
DEFAULT_RECIP_CELLS = 2
DEFAULT_CHARGE_PRECISION = 3


def compute_charges(
    structure,
    *,
    lambda_screening: float = DEFAULT_LAMBDA,
    h_electron_affinity: float = DEFAULT_H_ELECTRON_AFFINITY,
    eta: float = DEFAULT_ETA,
    real_cells: int = DEFAULT_REAL_CELLS,
    recip_cells: int = DEFAULT_RECIP_CELLS,
    charge_precision: int = DEFAULT_CHARGE_PRECISION,
    total_charge: float = 0.0,
) -> np.ndarray:
    """Return EQeq partial charges (one per site) for a pymatgen Structure.

    Faithful to EQeq's periodic ``ewald`` method with its default parameters.
    """
    from mofchecker_next._rust import eqeq_charges

    frac_coords = np.asarray(structure.frac_coords, dtype=float).tolist()
    a, b, c, alpha, beta, gamma = (float(v) for v in structure.lattice.parameters)

    electronegativity: list[float] = []
    hardness: list[float] = []
    for site in structure:
        params: AtomParameters = parameters_for(
            str(site.specie.symbol), h_electron_affinity
        )
        electronegativity.append(params.electronegativity)
        hardness.append(params.hardness)

    charges = eqeq_charges(
        frac_coords,
        [a, b, c],
        [alpha, beta, gamma],
        electronegativity,
        hardness,
        float(total_charge),
        float(lambda_screening),
        float(eta),
        int(real_cells),
        int(recip_cells),
        int(charge_precision),
    )
    return np.asarray(charges, dtype=float)


def high_charge_indices(
    structure,
    *,
    threshold: float = DEFAULT_CHARGE_THRESHOLD,
    **kwargs,
) -> list[int]:
    """Indices of sites whose absolute equilibrated charge exceeds ``threshold``."""
    charges = compute_charges(structure, **kwargs)
    return [i for i, q in enumerate(charges) if abs(q) > threshold]


def has_high_charges(
    structure,
    *,
    threshold: float = DEFAULT_CHARGE_THRESHOLD,
    **kwargs,
) -> bool:
    """True if any equilibrated charge exceeds ``threshold`` in magnitude.

    Mirrors MOFChecker 2.0's ``has_high_charges`` (which reports the inverse of an
    "are charges ok" check). The reference returns ``None`` when EQeq is missing;
    here the check always runs against the in-tree kernel.
    """
    return len(high_charge_indices(structure, threshold=threshold, **kwargs)) > 0


def check_high_charges(
    structure,
    *,
    threshold: float = DEFAULT_CHARGE_THRESHOLD,
    **kwargs,
) -> list[Diagnostic]:
    """Structured diagnostics for sites carrying unphysically high charges."""
    charges = compute_charges(structure, **kwargs)
    diagnostics: list[Diagnostic] = []
    for index, charge in enumerate(charges):
        if abs(charge) > threshold:
            diagnostics.append(
                Diagnostic(
                    check="high_charge",
                    severity="warning",
                    message="Equilibrated partial charge exceeds the plausible range",
                    atoms=[AtomRef(index=index, image=(0, 0, 0))],
                    values={
                        "charge": float(charge),
                        "threshold": float(threshold),
                    },
                )
            )
    return diagnostics
