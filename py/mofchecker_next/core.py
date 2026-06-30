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

VDW_H_RADIUS = 1.1
COVALENT_MEDIAN = 1.49


def _structure_from_file(path):
    """``Structure.from_file`` with pymatgen's benign CIF-rounding notice muted.

    pymatgen's CifParser warns whenever it snaps near-integer fractional
    coordinates to ideal values -- common for model-generated CIFs and harmless
    for the diagnostics. Scoped by message so genuine parse warnings still show.
    """
    import warnings

    from pymatgen.core import Structure

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r".*(fractional coordinates rounded to ideal values"
            r"|Issues encountered while parsing CIF).*",
        )
        return Structure.from_file(str(path))


def normalize_structure(obj):
    """Coerce a pymatgen Structure, ASE Atoms, or CIF path into a Structure."""
    from pymatgen.core import IStructure, Structure

    if isinstance(obj, (Structure, IStructure)):
        return obj
    if isinstance(obj, (str, Path)):
        return _structure_from_file(obj)
    if hasattr(obj, "get_chemical_symbols") and hasattr(obj, "get_positions"):
        from pymatgen.io.ase import AseAtomsAdaptor

        return AseAtomsAdaptor.get_structure(obj)
    raise TypeError(
        f"Unsupported structure input {type(obj)!r}; expected pymatgen Structure, "
        "ASE Atoms, or a CIF path."
    )

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
    ):
        """``distance_scale`` multiplies the bond-distance cutoffs used to build
        the neighbor graph (affects undercoordination, lone-molecule,
        connectivity, OMS, and the graph hashes). ``clash_scale`` multiplies the
        covalent-radius cutoffs used for atomic-overlap detection (affects
        ``has_atomic_overlaps``). Both default to ``1.0``, which reproduces
        MOFChecker exactly."""
        self.structure = structure
        self.metals = _co.METALS if metals is None else frozenset(str(m) for m in metals)
        self._method = method
        self._distance_scale = distance_scale
        self._clash_scale = clash_scale
        self._name = name
        self._path = path

    # -- constructors ------------------------------------------------------
    @classmethod
    def from_cif(cls, path, **kwargs):
        """Build from a CIF path (loaded with pymatgen ``Structure.from_file``)."""
        path = str(path)
        return cls(_structure_from_file(path), name=Path(path).stem, path=str(Path(path).resolve()), **kwargs)

    @classmethod
    def from_ase(cls, atoms, **kwargs):
        """Build from an ASE ``Atoms`` object."""
        from pymatgen.io.ase import AseAtomsAdaptor

        return cls(AseAtomsAdaptor.get_structure(atoms), **kwargs)

    @classmethod
    def from_structure(cls, structure, **kwargs):
        return cls(structure, **kwargs)

    # -- graph (built once) ------------------------------------------------
    @cached_property
    def graph(self):
        """The pymatgen StructureGraph (built once, reused by all checks)."""
        return _g.build_structure_graph(self.structure, self._method, distance_scale=self._distance_scale)

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
    def _overlap_indices(self) -> list[int]:
        from pymatgen.core import Element

        atomic_numbers = [int(site.specie.Z) for site in self.structure]
        radii_by_z = {
            int(Element(str(site.specie.symbol)).Z): _co.COVALENT_RADII.get(str(site.specie.symbol), COVALENT_MEDIAN)
            for site in self.structure
        }
        matrix = _geo.build_overlap_cutoff_matrix(atomic_numbers, radii_by_z, default_radius=COVALENT_MEDIAN)
        contacts = _geo.check_atomic_overlaps(self.structure, matrix, scale=self._clash_scale)
        idx = set()
        for d in contacts:
            for atom in d.atoms:
                idx.add(int(atom.index))
        return sorted(idx)

    def get_overlapping_indices(self) -> list[int]:
        return self._overlap_indices

    @property
    def has_atomic_overlaps(self) -> bool:
        return len(self._overlap_indices) > 0

    # -- coordination ------------------------------------------------------
    @cached_property
    def overvalent_c_indices(self) -> list[int]:
        return _g.overcoordinated_carbon_indices_from_structure(self.structure, self.metals, graph=self.graph)

    @cached_property
    def overcoordinated_n_indices(self) -> list[int]:
        return _g.overcoordinated_nitrogen_indices_from_structure(self.structure, self.metals, graph=self.graph)

    @cached_property
    def overvalent_h_indices(self) -> list[int]:
        return _geo.overcoordinated_hydrogen_indices(self.structure, VDW_H_RADIUS)

    @cached_property
    def undercoordinated_c_indices(self) -> list[int]:
        return _g.undercoordinated_carbon_indices_from_structure(
            self.structure, self.metals, _co.COVALENT_RADII, graph=self.graph
        )

    @cached_property
    def undercoordinated_n_indices(self) -> list[int]:
        return _g.undercoordinated_nitrogen_indices_from_structure(self.structure, self.metals, graph=self.graph)

    @cached_property
    def undercoordinated_rare_earth_indices(self) -> list[int]:
        return _g.undercoordinated_rare_earth_indices_from_structure(self.structure, graph=self.graph)

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
        return _g.floating_solvent_indices_from_structure(self.structure, graph=self.graph)

    @property
    def stray_atom_indices(self) -> list:
        """Detached finite components containing exactly one atom."""
        return [idx for idx in self.floating_solvent_indices if len(idx) == 1]

    @property
    def has_stray_atom(self) -> bool:
        return len(self.stray_atom_indices) > 0

    @property
    def lone_molecule_indices(self) -> list:
        """Detached finite components containing two or more atoms."""
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

    @property
    def has_high_charges(self) -> bool:
        from mofchecker_next.eqeq import has_high_charges

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
    def get_mof_descriptors(self, descriptors: Sequence[str] | None = None) -> "OrderedDict[str, object]":
        """Return an ordered dict of descriptor name -> value.

        Defaults to ``DEFAULT_DESCRIPTORS`` (metadata, symmetry, hashes, and all
        implemented diagnostics). Pass an explicit list to select a subset.
        """
        names = list(DEFAULT_DESCRIPTORS) if descriptors is None else list(descriptors)
        return OrderedDict((name, getattr(self, name)) for name in names)
