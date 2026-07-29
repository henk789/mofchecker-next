"""Drop-in ``MOFChecker``-compatible interface.

Exposes the same properties and ``get_mof_descriptors`` API as MOFChecker 2.0,
backed by the verified parity implementations in ``checks/`` and the Rust
kernels. The structure graph is built once and reused across all checks.

Not implemented (healing/correction and porosity are out of scope):
``adding_hydrogen``/``adding_linker`` raise ``NotImplementedError``; ``is_porous``
returns ``None`` (no bundled Zeo++), matching the reference when it cannot run.
"""

from __future__ import annotations

import base64
import hashlib
from collections import Counter, OrderedDict
from functools import cached_property
from pathlib import Path
from typing import Sequence

from mofchecker_next.checks import charge_oms as _co
from mofchecker_next.checks import composition as _comp
from mofchecker_next.checks import geometry as _geo
from mofchecker_next.checks import graph as _g
from mofchecker_next.profiles import DEFAULT_MODE, MODES, resolve_profile

VDW_H_RADIUS = 1.1
COVALENT_MEDIAN = 1.49
_UNSET = object()


class InputValidationError(ValueError):
    """Corrected-mode input violates the documented trust-boundary contract."""


def _structure_from_file(path, *, corrected: bool = False):
    """``Structure.from_file`` with pymatgen's benign CIF-rounding notice muted.

    pymatgen's CifParser warns whenever it snaps near-integer fractional
    coordinates to ideal values -- common for model-generated CIFs and harmless
    for the diagnostics. Scoped by message so genuine parse warnings still show.
    """
    import warnings

    if corrected:
        from pymatgen.io.cif import CifParser

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=r".*Issues encountered while parsing CIF.*")
            structures = CifParser(
                str(path), site_tolerance=0.0, frac_tolerance=0.0, check_cif=True,
            ).parse_structures(primitive=False, check_occu=True, on_error="raise")
        if len(structures) != 1:
            raise InputValidationError(f"expected one CIF structure, parsed {len(structures)}")
        return structures[0]

    from pymatgen.core import Structure

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r".*(fractional coordinates rounded to ideal values"
            r"|Issues encountered while parsing CIF).*",
        )
        return Structure.from_file(str(path))


def normalize_structure(obj, *, corrected: bool = False):
    """Coerce a pymatgen Structure, ASE Atoms, or CIF path into a Structure."""
    from pymatgen.core import IStructure, Structure

    if isinstance(obj, (Structure, IStructure)):
        return obj
    if isinstance(obj, (str, Path)):
        return _structure_from_file(obj, corrected=corrected)
    if hasattr(obj, "get_chemical_symbols") and hasattr(obj, "get_positions"):
        from pymatgen.io.ase import AseAtomsAdaptor

        return AseAtomsAdaptor.get_structure(obj)
    raise TypeError(
        f"Unsupported structure input {type(obj)!r}; expected pymatgen Structure, "
        "ASE Atoms, or a CIF path."
    )


def validate_corrected_structure(structure) -> None:
    """Reject inputs for which corrected diagnostics would fabricate certainty."""
    import numpy as np

    lattice = np.asarray(structure.lattice.matrix, dtype=float)
    frac = np.asarray(structure.frac_coords, dtype=float)
    cart = np.asarray(structure.cart_coords, dtype=float)
    if len(structure) == 0:
        raise InputValidationError("structure contains no sites")
    if lattice.shape != (3, 3) or not np.isfinite(lattice).all():
        raise InputValidationError("lattice contains non-finite values")
    determinant = float(np.linalg.det(lattice))
    if not np.isfinite(determinant) or abs(determinant) <= 1e-8:
        raise InputValidationError("lattice is singular or has negligible volume")
    if not np.isfinite(frac).all() or not np.isfinite(cart).all():
        raise InputValidationError("site coordinates contain NaN or infinity")
    for index, site in enumerate(structure):
        if not site.is_ordered or len(site.species) != 1 or abs(float(site.species.num_atoms) - 1.0) > 1e-8:
            raise InputValidationError(f"site {index} has disorder or partial occupancy")
        try:
            atomic_number = int(site.specie.Z)
        except Exception as exc:
            raise InputValidationError(f"site {index} is not a supported chemical element") from exc
        if not 1 <= atomic_number <= 118:
            raise InputValidationError(f"site {index} has invalid atomic number {atomic_number}")


# Descriptors returned by get_mof_descriptors() by default: metadata, symmetry,
# graph hashes, and every implemented diagnostic. Excludes the healing
# descriptors (adding_hydrogen/adding_linker) which raise NotImplementedError.
DEFAULT_DESCRIPTORS = (
    "name", "formula", "density", "spacegroup_symbol", "spacegroup_number",
    "graph_hash", "undecorated_graph_hash", "scaffold_hash",
    "undecorated_scaffold_hash", "symmetry_hash",
    "has_carbon", "has_hydrogen", "has_nitrogen", "has_metal", "metal_number",
    "has_atomic_overlaps", "has_overcoordinated_c", "has_overcoordinated_n",
    "has_overcoordinated_h", "has_undercoordinated_c", "has_undercoordinated_n",
    "has_undercoordinated_rare_earth", "has_undercoordinated_alkali_alkaline",
    "has_stray_atom", "has_lone_molecule", "has_3d_connected_graph", "has_suspicious_terminal_oxo",
    "has_geometrically_exposed_metal", "possible_charged_fused_ring",
    "positive_charge_from_linkers", "negative_charge_from_linkers",
    "has_high_charges", "has_oms", "is_porous",
)

NEXT_DEFAULT_DESCRIPTORS = DEFAULT_DESCRIPTORS + (
    "max_abs_eqeq_charge", "eqeq_charge_sum", "eqeq_expected_total_charge",
    "eqeq_charge_residual", "eqeq_charge_threshold",
)

LEGACY_DEFAULT_DESCRIPTORS = (
    "name", "graph_hash", "undecorated_graph_hash", "decorated_scaffold_hash",
    "undecorated_scaffold_hash", "symmetry_hash", "formula", "path", "density",
    "has_carbon", "has_hydrogen", "has_atomic_overlaps", "has_overcoordinated_c",
    "has_overcoordinated_n", "has_overcoordinated_h", "has_undercoordinated_c",
    "has_undercoordinated_n", "has_undercoordinated_rare_earth", "has_metal",
    "has_lone_molecule", "has_high_charges", "is_porous",
    "has_suspicicious_terminal_oxo", "has_undercoordinated_alkali_alkaline",
    "has_geometrically_exposed_metal", "has_3d_connected_graph",
)


class MOFChecker:
    """MOFChecker-compatible diagnostics for a single structure."""

    def __init__(
        self,
        structure,
        *,
        metals=None,
        method: str = "vesta",
        distance_scale: float = 1.0,
        clash_scale: float = 1.0,
        name=None,
        path=None,
        mode: str = DEFAULT_MODE,
        symprec=_UNSET,
        angle_tolerance=_UNSET,
        primitive=_UNSET,
    ):
        """``distance_scale`` multiplies the bond-distance cutoffs used to build
        the neighbor graph (affects undercoordination, lone-molecule,
        connectivity, OMS, and the graph hashes). ``clash_scale`` multiplies the
        covalent-radius cutoffs used for atomic-overlap detection (affects
        ``has_atomic_overlaps``). Both default to ``1.0``, which reproduces
        MOFChecker exactly."""
        profile = resolve_profile(mode)
        if profile.legacy_input:
            from pymatgen.core import Structure
            from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

            symprec = 0.5 if symprec is _UNSET else symprec
            angle_tolerance = 5 if angle_tolerance is _UNSET else angle_tolerance
            primitive = True if primitive is _UNSET else primitive
            if symprec is not None or angle_tolerance is not None:
                try:
                    structure = SpacegroupAnalyzer(
                        structure, symprec=symprec, angle_tolerance=angle_tolerance
                    ).get_symmetrized_structure()
                except TypeError:
                    pass
            if primitive:
                structure = structure.get_primitive_structure()
            if isinstance(structure, Structure):
                from pymatgen.core import IStructure

                structure = IStructure.from_sites(structure)

        if profile.corrected_input:
            validate_corrected_structure(structure)
        if not isinstance(method, str) or not method:
            raise ValueError("method must be a non-empty string")
        for label, value in (("distance_scale", distance_scale), ("clash_scale", clash_scale)):
            if not isinstance(value, (int, float)) or not 0 < float(value) < float("inf"):
                raise ValueError(f"{label} must be finite and positive")

        self.structure = structure
        self.mode = mode
        self.profile = profile
        self.metals = _co.METALS if metals is None else frozenset(str(m) for m in metals)
        self._method = method
        self._distance_scale = distance_scale
        self._clash_scale = clash_scale
        self._name = name
        self._path = path

    # -- constructors ------------------------------------------------------
    @classmethod
    def from_cif(cls, path, **kwargs):
        """Build from a CIF path using the selected reference mode's defaults."""
        path = str(path)
        if resolve_profile(kwargs.get("mode", DEFAULT_MODE)).legacy_input:
            import warnings
            from pymatgen.io.cif import CifParser

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                structure = CifParser(path).get_structures()[0]
            kwargs.setdefault("symprec", 0.5)
            kwargs.setdefault("angle_tolerance", 5)
            kwargs.setdefault("primitive", False)
        else:
            structure = _structure_from_file(path, corrected=resolve_profile(kwargs.get("mode", DEFAULT_MODE)).corrected_input)
        return cls(structure, name=Path(path).stem, path=str(Path(path).resolve()), **kwargs)

    @classmethod
    def from_ase(cls, atoms, **kwargs):
        """Build from an ASE ``Atoms`` object."""
        from pymatgen.io.ase import AseAtomsAdaptor

        if resolve_profile(kwargs.get("mode", DEFAULT_MODE)).legacy_input:
            kwargs.setdefault("symprec", 0.5)
            kwargs.setdefault("angle_tolerance", 5)
            kwargs.setdefault("primitive", False)
        return cls(AseAtomsAdaptor.get_structure(atoms), **kwargs)

    @classmethod
    def from_structure(cls, structure, **kwargs):
        return cls(structure, **kwargs)

    # -- graph (built once) ------------------------------------------------
    @cached_property
    def graph(self):
        """The pymatgen StructureGraph (built once, reused by all checks)."""
        structure = self.structure
        if self.profile.corrected_input:
            from pymatgen.core import Structure

            # Canonical computational view: preserve cell/index order while making
            # graph construction invariant to integer translations of sites.
            structure = Structure.from_sites(structure, to_unit_cell=True)
        return _g.build_structure_graph(structure, self._method, distance_scale=self._distance_scale)

    # -- metadata ----------------------------------------------------------
    @property
    def name(self) -> str | None:
        return self._name

    @property
    def path(self) -> str | None:
        return self._path

    @property
    def formula(self) -> str:
        return self.structure.composition.formula

    @property
    def density(self) -> float:
        return float(self.structure.density)

    @property
    def volume(self) -> float:
        return float(self.structure.volume)

    # -- symmetry ----------------------------------------------------------
    @cached_property
    def _symmetrized(self):
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        return SpacegroupAnalyzer(self.structure).get_symmetrized_structure()

    @property
    def spacegroup_symbol(self) -> str:
        return self._symmetrized.spacegroup.int_symbol

    @property
    def spacegroup_number(self) -> int:
        return int(self._symmetrized.spacegroup.int_number)

    @cached_property
    def symmetry_hash(self) -> str:
        """Hash of the symmetrized structure (Wyckoff-letter set + spacegroup).

        Note: this is *deterministic* (Wyckoff letters are sorted before
        hashing). The reference MOFChecker hashes ``tuple(set(...))`` without
        sorting, so its value depends on Python's per-process string-hash
        randomization and is not reproducible across runs; our values will not
        match it but are stable.
        """
        sym = self._symmetrized
        wyckoff = tuple(sorted(set(sym.wyckoff_letters)))
        hasher = hashlib.sha256()
        hasher.update(repr(wyckoff).encode())
        return base64.b64encode(hasher.digest()).decode() + str(sym.spacegroup.int_number)

    # -- graph hashes ------------------------------------------------------
    @cached_property
    def graph_hash(self) -> str:
        from structuregraph_helpers.hash import decorated_graph_hash

        return decorated_graph_hash(self.graph, lqg=False)

    @cached_property
    def undecorated_graph_hash(self) -> str:
        from structuregraph_helpers.hash import undecorated_graph_hash

        return undecorated_graph_hash(self.graph, lqg=False)

    @cached_property
    def scaffold_hash(self) -> str:
        from structuregraph_helpers.hash import decorated_scaffold_hash

        return decorated_scaffold_hash(self.graph, lqg=False)

    @cached_property
    def undecorated_scaffold_hash(self) -> str:
        from structuregraph_helpers.hash import undecorated_scaffold_hash

        return undecorated_scaffold_hash(self.graph, lqg=False)

    # -- composition -------------------------------------------------------
    @property
    def has_carbon(self) -> bool:
        return _comp.has_element(self.structure, "C")

    @property
    def has_hydrogen(self) -> bool:
        return _comp.has_element(self.structure, "H")

    @property
    def has_nitrogen(self) -> bool:
        return _comp.has_element(self.structure, "N")

    @property
    def has_metal(self) -> bool:
        return _comp.has_metal(self.structure, self.metals)

    @property
    def metal_number(self) -> int:
        return _comp.metal_number(self.structure, self.metals)

    # -- atomic overlaps ---------------------------------------------------
    @cached_property
    def _overlap_diagnostics(self):
        from pymatgen.core import Element

        atomic_numbers = [int(site.specie.Z) for site in self.structure]
        radii_by_z = {
            int(Element(str(site.specie.symbol)).Z): _co.COVALENT_RADII.get(str(site.specie.symbol), COVALENT_MEDIAN)
            for site in self.structure
        }
        matrix = _geo.build_overlap_cutoff_matrix(atomic_numbers, radii_by_z, default_radius=COVALENT_MEDIAN)
        return _geo.check_atomic_overlaps(self.structure, matrix, scale=self._clash_scale)

    @cached_property
    def _overlap_indices(self) -> list[int]:
        idx = set()
        for d in self._overlap_diagnostics:
            for atom in d.atoms:
                idx.add(int(atom.index))
        return sorted(idx)

    @property
    def overlap_diagnostics(self):
        return self._overlap_diagnostics

    def get_overlapping_indices(self) -> list[int]:
        return self._overlap_indices

    @property
    def has_atomic_overlaps(self) -> bool:
        return len(self._overlap_indices) > 0

    # -- coordination ------------------------------------------------------
    @cached_property
    def overvalent_c_indices(self) -> list[int]:
        return _g.overcoordinated_carbon_indices_from_structure(
            self.structure, self.metals, graph=self.graph, exclude_boron=self.mode != "0.9.6"
        )

    @cached_property
    def overcoordinated_n_indices(self) -> list[int]:
        return _g.overcoordinated_nitrogen_indices_from_structure(self.structure, self.metals, graph=self.graph)

    @cached_property
    def overvalent_h_indices(self) -> list[int]:
        return _geo.overcoordinated_hydrogen_indices(self.structure, VDW_H_RADIUS)

    @cached_property
    def undercoordinated_c_indices(self) -> list[int]:
        if self.mode == "0.9.6":
            return _g.undercoordinated_carbon_indices_v096_from_structure(
                self.structure, graph=self.graph
            )
        return _g.undercoordinated_carbon_indices_from_structure(
            self.structure, self.metals, _co.COVALENT_RADII, graph=self.graph,
            use_connected_site_vectors=self.profile.corrected_carbon_images,
        )

    @cached_property
    def undercoordinated_n_indices(self) -> list[int]:
        if self.mode == "0.9.6":
            return _g.undercoordinated_nitrogen_indices_v096_from_structure(
                self.structure, graph=self.graph
            )
        return _g.undercoordinated_nitrogen_indices_from_structure(self.structure, self.metals, graph=self.graph)

    @cached_property
    def undercoordinated_rare_earth_indices(self) -> list[int]:
        # 0.9.6 used pymatgen's is_rare_earth_metal (no Sc/Y); 2.0 switched to
        # is_rare_earth, which also counts Sc and Y.
        return _g.undercoordinated_rare_earth_indices_from_structure(
            self.structure, graph=self.graph, include_sc_y=self.mode != "0.9.6"
        )

    @cached_property
    def _undercoordinated_alkali_alkaline_indices(self) -> list[int]:
        return _g.undercoordinated_alkali_alkaline_indices_from_structure(self.structure, graph=self.graph)

    @property
    def has_overcoordinated_c(self) -> bool:
        return len(self.overvalent_c_indices) > 0

    @property
    def has_overcoordinated_n(self) -> bool:
        return len(self.overcoordinated_n_indices) > 0

    @property
    def has_overcoordinated_h(self) -> bool:
        return len(self.overvalent_h_indices) > 0

    @property
    def has_undercoordinated_c(self) -> bool:
        return len(self.undercoordinated_c_indices) > 0

    @property
    def has_undercoordinated_n(self) -> bool:
        return len(self.undercoordinated_n_indices) > 0

    @property
    def has_undercoordinated_rare_earth(self) -> bool:
        return len(self.undercoordinated_rare_earth_indices) > 0

    @property
    def has_undercoordinated_alkali_alkaline(self) -> bool:
        return len(self._undercoordinated_alkali_alkaline_indices) > 0

    # -- floating solvent / connectivity ----------------------------------
    @cached_property
    def floating_solvent_indices(self) -> list:
        """All finite detached components (old lone_molecule_indices behavior)."""
        if self.mode == "0.9.6":
            return _g.floating_solvent_indices_v096_from_structure(self.structure, graph=self.graph)
        return _g.floating_solvent_indices_from_structure(self.structure, graph=self.graph)

    @property
    def stray_atom_indices(self) -> list:
        """Detached finite components containing exactly one atom."""
        if self.mode == "0.9.6":
            return []
        return [idx for idx in self.floating_solvent_indices if len(idx) == 1]

    @property
    def has_stray_atom(self) -> bool:
        return len(self.stray_atom_indices) > 0

    @property
    def lone_molecule_indices(self) -> list:
        """Detached finite components; 0.9.6 includes single stray atoms here."""
        if self.mode == "0.9.6":
            return self.floating_solvent_indices
        return [idx for idx in self.floating_solvent_indices if len(idx) >= 2]

    @property
    def has_lone_molecule(self) -> bool:
        return len(self.lone_molecule_indices) > 0

    @property
    def has_3d_connected_graph(self) -> bool:
        return _g.is_3d_connected_graph_from_structure(self.structure, graph=self.graph)

    # -- metal-site checks -------------------------------------------------
    @cached_property
    def suspicious_terminal_oxo_indices(self) -> list[int]:
        return _g.false_oxo_indices_from_structure(self.structure, self.metals, graph=self.graph)

    @property
    def has_suspicious_terminal_oxo(self) -> bool:
        return len(self.suspicious_terminal_oxo_indices) > 0

    @cached_property
    def geometrically_exposed_metal_indice(self) -> list[int]:
        return _g.geometrically_exposed_metal_indices_from_structure(self.structure, self.metals, graph=self.graph)

    @property
    def has_geometrically_exposed_metal(self) -> bool:
        return len(self.geometrically_exposed_metal_indice) > 0

    @cached_property
    def oms_indice(self) -> list[int]:
        return _co.oms_indices(self.structure, self.graph)

    @property
    def has_oms(self) -> bool:
        return _co.has_oms(self.structure, self.graph)

    # -- charge ------------------------------------------------------------
    @cached_property
    def _clean_cycles(self):
        # ponytail: NetworkX cycle enumeration is exact but expensive; compute it once per structure.
        return _co._clean_cycles(self.graph)

    @property
    def possible_charged_fused_ring(self) -> bool:
        return len(_co.fused_ring_indices(self.structure, self.graph, cycles=self._clean_cycles)) > 0

    @property
    def positive_charge_from_linkers(self) -> int:
        return len(_co.positive_charge_indices(self.structure, self.graph, cycles=self._clean_cycles))

    @property
    def negative_charge_from_linkers(self) -> int:
        return len(_co.negative_charge_indices(self.structure, self.graph, cycles=self._clean_cycles))

    @cached_property
    def eqeq_charges(self) -> tuple[float, ...]:
        """JSON-safe true-element Rust EQeq charges for the corrected profile."""
        if self.profile.total_charge is None:
            raise AttributeError("this profile does not declare an EQeq total charge")
        from mofchecker_next.eqeq import compute_charges

        charges = tuple(float(q) for q in compute_charges(
            self.structure,
            total_charge=self.profile.total_charge,
            reference_cif_labels=False,
        ))
        residual = abs(sum(charges) - self.profile.total_charge)
        if residual > self.profile.eqeq_charge_sum_tolerance:
            raise ValueError(
                f"EQeq charge sum {sum(charges):.12g} differs from requested "
                f"{self.profile.total_charge:.12g} by more than "
                f"{self.profile.eqeq_charge_sum_tolerance:g}"
            )
        return charges

    @property
    def max_abs_eqeq_charge(self) -> float:
        return max((abs(q) for q in self.eqeq_charges), default=0.0)

    @property
    def eqeq_charge_sum(self) -> float:
        return sum(self.eqeq_charges)

    @property
    def eqeq_expected_total_charge(self) -> float:
        if self.profile.total_charge is None:
            raise AttributeError("this profile does not declare an EQeq total charge")
        return self.profile.total_charge

    @property
    def eqeq_charge_residual(self) -> float:
        return abs(sum(self.eqeq_charges) - self.eqeq_expected_total_charge)

    @property
    def eqeq_charge_threshold(self) -> float:
        return self.profile.eqeq_threshold

    @property
    def has_high_charges(self) -> bool:
        from mofchecker_next.eqeq import has_high_charges

        if self.mode == "0.9.6":
            # 0.9.6 threshold, and EQeq fed through MOFChecker's CIF round-trip.
            return has_high_charges(self.structure, threshold=3.0, reference_cif_labels=True)
        if self.profile.corrected_input:
            return self.max_abs_eqeq_charge > self.eqeq_charge_threshold
        return has_high_charges(self.structure)

    # -- out of scope ------------------------------------------------------
    @property
    def is_porous(self):
        """Porosity (Zeo++) is not bundled; returns None like the reference when
        it cannot run."""
        return None

    @property
    def adding_hydrogen(self):
        raise NotImplementedError(
            "Hydrogen addition (healing) is out of scope for mofchecker-next; "
            "use the original MOFChecker for correction workflows."
        )

    @property
    def adding_linker(self):
        raise NotImplementedError(
            "Linker addition (healing) is out of scope for mofchecker-next."
        )

    # -- descriptor dict ---------------------------------------------------
    @property
    def decorated_scaffold_hash(self) -> str:
        """0.9.6 name for ``scaffold_hash``."""
        return self.scaffold_hash

    @property
    def has_suspicicious_terminal_oxo(self) -> bool:
        """0.9.6's misspelled public property."""
        return self.has_suspicious_terminal_oxo

    def get_mof_descriptors(self, descriptors: Sequence[str] | None = None) -> "OrderedDict[str, object]":
        """Return an ordered dict of descriptor name -> value.

        Defaults to ``DEFAULT_DESCRIPTORS`` (metadata, symmetry, hashes, and all
        implemented diagnostics). Pass an explicit list to select a subset.
        """
        defaults = {
            "legacy": LEGACY_DEFAULT_DESCRIPTORS,
            "modern": DEFAULT_DESCRIPTORS,
            "next": NEXT_DEFAULT_DESCRIPTORS,
        }[self.profile.descriptor_set]
        names = list(defaults) if descriptors is None else list(descriptors)
        return OrderedDict((name, getattr(self, name)) for name in names)
