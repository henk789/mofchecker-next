from __future__ import annotations

from dataclasses import asdict
from typing import Mapping, Sequence

import numpy as np

from mofchecker_next.diagnostics import AtomRef, Diagnostic


def minimum_image_distance_py(frac_i, frac_j, lattice):
    """Return minimum-image displacement for fractional coordinates.

    ``delta_frac`` is r_ij in fractional coordinates: it points from atom i to
    the selected periodic image of atom j.

    The starting image follows ``image = -np.round(frac_j - frac_i)``.
    NumPy rounds exact half-integers to the nearest even integer, so a
    displacement of exactly ``+0.5`` or ``-0.5`` remains in the central image.

    For skewed non-orthogonal cells, component-wise fractional wrapping can miss
    the shortest Cartesian image. To avoid that MIC pitfall, this checks the 26
    neighboring images around the rounded image and keeps the first exact tie.
    """
    frac_i = np.asarray(frac_i, dtype=float)
    frac_j = np.asarray(frac_j, dtype=float)
    lattice = np.asarray(lattice, dtype=float)

    d = frac_j - frac_i
    base_image = -np.round(d).astype(int)

    best_image = base_image.copy()
    best_delta_frac = d + best_image
    best_delta_cart = best_delta_frac @ lattice
    best_distance = float(np.linalg.norm(best_delta_cart))

    offsets = [(0, 0, 0)]
    offsets.extend(
        (x, y, z)
        for x in (-1, 0, 1)
        for y in (-1, 0, 1)
        for z in (-1, 0, 1)
        if (x, y, z) != (0, 0, 0)
    )
    for offset in offsets[1:]:
        image = base_image + np.array(offset, dtype=int)
        delta_frac = d + image
        delta_cart = delta_frac @ lattice
        distance = float(np.linalg.norm(delta_cart))
        if distance < best_distance:
            best_image = image
            best_delta_frac = delta_frac
            best_distance = distance

    return {
        "delta_frac": best_delta_frac.tolist(),
        "image": best_image.tolist(),
        "distance": best_distance,
    }


def find_short_contacts_py(
    frac_coords,
    atomic_numbers,
    lattice,
    cutoff_matrix,
    scale: float = 1.0,
):
    """Brute-force short-contact detection using supplied cutoffs."""
    frac_coords = np.asarray(frac_coords, dtype=float)
    atomic_numbers = [int(number) for number in atomic_numbers]
    lattice = np.asarray(lattice, dtype=float)

    contacts = []
    for i in range(len(frac_coords) - 1):
        zi = atomic_numbers[i]
        for j in range(i + 1, len(frac_coords)):
            zj = atomic_numbers[j]
            cutoff = float(cutoff_matrix[zi][zj]) * scale
            image_distance = minimum_image_distance_py(frac_coords[i], frac_coords[j], lattice)
            distance = image_distance["distance"]
            if distance < cutoff:
                contacts.append(
                    {
                        "i": i,
                        "j": j,
                        "image_j": [int(v) for v in image_distance["image"]],
                        "distance": distance,
                        "cutoff": cutoff,
                    }
                )

    contacts.sort(key=lambda c: (c["i"], c["j"], c["image_j"], c["distance"]))
    return contacts


def find_neighbor_candidates_py(frac_coords, lattice, cutoff: float):
    """Brute-force scalar-cutoff neighbor candidates using the same MIC."""
    frac_coords = np.asarray(frac_coords, dtype=float)
    lattice = np.asarray(lattice, dtype=float)
    candidates = []
    for i in range(len(frac_coords) - 1):
        for j in range(i + 1, len(frac_coords)):
            image_distance = minimum_image_distance_py(frac_coords[i], frac_coords[j], lattice)
            distance = image_distance["distance"]
            if distance < cutoff:
                candidates.append(
                    {
                        "i": i,
                        "j": j,
                        "image_j": [int(v) for v in image_distance["image"]],
                        "distance": distance,
                    }
                )
    candidates.sort(key=lambda c: (c["i"], c["j"], c["image_j"], c["distance"]))
    return candidates


def build_overlap_cutoff_matrix(
    atomic_numbers: Sequence[int],
    radii_by_atomic_number: Mapping[int, float],
    default_radius: float | None = None,
) -> list[list[float]]:
    """Build an overlap cutoff matrix from caller-supplied covalent radii.

    MOFChecker 2.0 uses ``min(covalent_radius_i, covalent_radius_j)`` for its
    atomic-overlap criterion. This helper implements only that matrix operation;
    it does not infer radii and keeps chemistry data on the Python side.
    """
    atomic_numbers = [int(number) for number in atomic_numbers]
    if not atomic_numbers:
        return []

    size = max(atomic_numbers) + 1
    matrix = [[0.0 for _ in range(size)] for _ in range(size)]
    radii = {}
    for atomic_number in set(atomic_numbers):
        if atomic_number in radii_by_atomic_number:
            radii[atomic_number] = float(radii_by_atomic_number[atomic_number])
        elif default_radius is not None:
            radii[atomic_number] = float(default_radius)
        else:
            raise KeyError(f"Missing covalent radius for atomic number {atomic_number}")

    for zi, ri in radii.items():
        for zj, rj in radii.items():
            matrix[zi][zj] = min(ri, rj)
    return matrix


def structure_to_arrays(structure):
    """Extract array-like geometry inputs from a pymatgen Structure."""
    frac_coords = np.asarray(structure.frac_coords, dtype=float).tolist()
    atomic_numbers = [int(site.specie.Z) for site in structure]
    lattice = np.asarray(structure.lattice.matrix, dtype=float).tolist()
    return frac_coords, atomic_numbers, lattice


def check_atomic_overlaps(structure, cutoff_matrix, scale: float = 1.0) -> list[Diagnostic]:
    """Run the Rust short-contact kernel and return structured diagnostics."""
    from mofchecker_next._rust import find_short_contacts

    frac_coords, atomic_numbers, lattice = structure_to_arrays(structure)
    contacts = find_short_contacts(frac_coords, atomic_numbers, lattice, cutoff_matrix, scale)
    diagnostics = []
    for contact in contacts:
        diagnostics.append(
            Diagnostic(
                check="atomic_overlap",
                severity="error",
                message="Unphysically short interatomic distance",
                atoms=[
                    AtomRef(index=int(contact["i"]), image=(0, 0, 0)),
                    AtomRef(index=int(contact["j"]), image=tuple(int(v) for v in contact["image_j"])),
                ],
                values={
                    "distance_angstrom": float(contact["distance"]),
                    "cutoff_angstrom": float(contact["cutoff"]),
                },
            )
        )
    return diagnostics


def diagnostics_to_dicts(diagnostics: Sequence[Diagnostic]) -> list[dict]:
    return [asdict(diagnostic) for diagnostic in diagnostics]


def overcoordinated_hydrogen_indices(structure, vdw_h_radius: float) -> list[int]:
    """Match MOFChecker 2.0's H overcoordination neighbor-count rule.

    MOFChecker 2.0 flags hydrogen atoms with more than one pymatgen neighbor
    within the van-der-Waals radius of H. The radius is supplied by the caller.
    """
    flagged = []
    for index, site in enumerate(structure):
        if int(site.specie.Z) != 1:
            continue
        if len(structure.get_neighbors(site, float(vdw_h_radius))) > 1:
            flagged.append(index)
    return flagged


def check_overcoordinated_hydrogen(structure, vdw_h_radius: float) -> list[Diagnostic]:
    diagnostics = []
    for index in overcoordinated_hydrogen_indices(structure, vdw_h_radius):
        diagnostics.append(
            Diagnostic(
                check="overcoordinated_hydrogen",
                severity="error",
                message="Hydrogen has more than one neighbor within the H van-der-Waals radius",
                atoms=[AtomRef(index=index)],
                values={"neighbor_cutoff_angstrom": float(vdw_h_radius)},
            )
        )
    return diagnostics
