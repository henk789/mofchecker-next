"""Linker-charge, fused-ring, and open-metal-site diagnostics.

Faithful functional ports of MOFChecker 2.0's graph-based heuristics
(`positive_charge.py`, `negative_charge.py`, `fused_ring.py`, `oms/`). They run
on a pymatgen ``StructureGraph`` built the same way the reference builds it
(`structuregraph_helpers`, "vesta" method), and reuse the same library calls
(`construct_clean_graph`, `get_cn`, `nx.simple_cycles`, `LocalStructOrderParams`)
so the results match the reference. The metal set and covalent radii mirror the
reference tables (METALS; Cordero et al. 2008 covalent radii).
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from pymatgen.util.coord import get_angle

from mofchecker_next.checks.graph import build_structure_graph

# Metals per MOFChecker 2.0's definitions.METALS (transition metals, lanthanides,
# actinides, alkali/alkaline-earth, and Al, Ga, In, Tl, Ge, Sn, Sb, Bi, Po, ...).
METALS = frozenset(
    "Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge Rb Sr Y Zr Nb Mo Tc "
    "Ru Rh Pd Ag Cd In Sn Sb Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
    "Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es "
    "Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn".split()
)

_HALOGENS = frozenset({"F", "Cl", "Br", "I"})

# Cordero et al. 2008 covalent radii (DOI:10.1039/B801115J), plain element symbols.
COVALENT_RADII = {
    "H": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76,
    "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58, "Na": 1.66, "Mg": 1.41,
    "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Ar": 1.06,
    "K": 2.03, "Ca": 1.76, "Sc": 1.7, "Ti": 1.6, "V": 1.53, "Cr": 1.39,
    "Mn": 1.61, "Fe": 1.52, "Co": 1.5, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
    "Ga": 1.22, "Ge": 1.2, "As": 1.19, "Se": 1.2, "Br": 1.2, "Kr": 1.16,
    "Rb": 2.2, "Sr": 1.95, "Y": 1.9, "Zr": 1.5, "Nb": 1.64, "Mo": 1.54,
    "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44,
    "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "Xe": 1.4,
    "Cs": 2.44, "Ba": 2.15, "La": 2.07, "Ce": 2.04, "Pr": 2.03, "Nd": 2.01,
    "Pm": 1.99, "Sm": 1.98, "Eu": 1.98, "Gd": 1.96, "Tb": 1.94, "Dy": 1.92,
    "Ho": 1.92, "Er": 1.89, "Tm": 1.9, "Yb": 1.87, "Lu": 1.87, "Hf": 1.75,
    "Ta": 1.7, "W": 1.62, "Re": 1.51, "Os": 1.44, "Ir": 1.41, "Pt": 1.36,
    "Au": 1.36, "Hg": 1.32, "Tl": 1.45, "Pb": 1.46, "Bi": 1.48, "Po": 1.4,
    "At": 1.5, "Rn": 1.5, "Fr": 2.6, "Ra": 2.21, "Ac": 2.15, "Th": 2.06,
    "Pa": 2.0, "U": 1.96, "Np": 1.9, "Pu": 1.87, "Am": 1.8, "Cm": 1.69,
}

# Order-parameter thresholds for OMS detection (MOFChecker 2.0 oms/definitions.py).
OP_DEF = {
    4: {"names": ["sq_plan", "sq", "see_saw_rect", "tet", "tri_pyr"],
        "weights": [0.2, 0.1, 0.1, 0.5, 0.5], "open": [0, 1, 2, 4]},
    5: {"names": ["pent_plan", "sq_pyr", "tri_bipyr"], "weights": [1, 0.5, 0.5], "open": [0, 1]},
    6: {"names": ["pent_pyr", "oct"], "weights": [0.3, 0.7], "open": [0]},
    7: {"names": ["hex_pyr", "pent_bipyr"], "weights": [0.7, 0.3], "open": [0]},
    8: {"names": ["hex_bipyr"], "weights": [1], "open": None},
}


# ---------------------------------------------------------------------------
# Small helpers mirroring the reference get_indices utilities.
# ---------------------------------------------------------------------------
def _is_metal(site) -> bool:
    return str(site.specie) in METALS


def _species_list(structure) -> tuple[str, ...]:
    species = getattr(structure, "_mofchecker_next_species", None)
    if species is None:
        species = structure._mofchecker_next_species = tuple(str(site.specie) for site in structure)
    return species


def _is_halogen(site) -> bool:
    return str(site.specie) in _HALOGENS


def _num_neighbor_metal(neighbors) -> int:
    return sum(1 for n in neighbors if _is_metal(n.site))


def _num_neighbor_halogen(neighbors) -> int:
    return sum(1 for n in neighbors if _is_halogen(n.site))


def _non_metal_neighbor(neighbors):
    for n in neighbors:
        if not _is_metal(n.site):
            return n
    return None


def _non_metal_neighbors(neighbors):
    return [n for n in neighbors if not _is_metal(n.site)]


def _connected(graph, index: int):
    cache = getattr(graph, "_mofchecker_next_connected_sites", None)
    if cache is None:
        cache = graph._mofchecker_next_connected_sites = {}
    if index not in cache:
        cache[index] = graph.get_connected_sites(index)
    return cache[index]


def _cn(graph, index: int) -> int:
    return len(_connected(graph, index))


def _metal_neighbor_count(graph, index: int) -> int:
    cache = getattr(graph, "_mofchecker_next_metal_neighbor_count", None)
    if cache is None:
        cache = graph._mofchecker_next_metal_neighbor_count = {}
    if index not in cache:
        cache[index] = _num_neighbor_metal(_connected(graph, index))
    return cache[index]


def _species(structure, index: int) -> str:
    return _species_list(structure)[index]


def _indices(structure, symbol: str) -> list[int]:
    return [i for i, site in enumerate(structure) if str(site.specie) == symbol]


def _clean_cycles(graph) -> list[list[int]]:
    from mofchecker_next._rust import bounded_simple_cycles_undirected

    # Same simple graph as structuregraph_helpers.construct_clean_graph(...),
    # but skip NetworkX object construction and enumerate bounded cycles in Rust.
    edges = {(min(u, v), max(u, v)) for u, v in graph.graph.edges() if u != v}
    return bounded_simple_cycles_undirected(len(graph.structure), sorted(edges), 16)


def _nonmetal_cycles_by_node(structure, cycles, length: int) -> list[list[list[int]]]:
    by_node: list[list[list[int]]] = [[] for _ in structure]
    is_nonmetal = [sp not in METALS for sp in _species_list(structure)]
    for cycle in cycles:
        if len(cycle) == length and all(is_nonmetal[i] for i in cycle):
            for i in cycle:
                by_node[i].append(cycle)
    return by_node


def get_angle_between_site_and_neighbors(site, neighbors) -> float:
    """Minimum angle between a site and its two neighbors (degrees)."""
    vec_1 = site.coords - neighbors[1].site.coords
    vec_2 = site.coords - neighbors[0].site.coords
    return get_angle(vec_1, vec_2)


def _guess_underbound_nitrogen_cn2(
    structure,
    site_index: int,
    neighbors: list,
    connected_sites_a: list,
    connected_sites_b: list,
    tolerance: float = 25,
) -> bool:
    """Port of MOFChecker 2.0's CN-2 nitrogen under-coordination heuristic."""

    def vector_angle_to_plane(a, b, c, d):
        a, b, c, d = map(np.array, [a, b, c, d])
        ab = b - a
        ac = c - a
        normal_vector = np.cross(ab, ac)
        dc = c - d
        dot_product = np.dot(dc, normal_vector)
        magnitude_dc = np.linalg.norm(dc)
        magnitude_normal = np.linalg.norm(normal_vector)
        angle_rad = np.arcsin(dot_product / (magnitude_dc * magnitude_normal))
        return abs(np.degrees(angle_rad))

    angle = get_angle_between_site_and_neighbors(structure[site_index], neighbors)
    num_h = 0
    bond_lengths = np.array(
        [
            structure.get_distance(site_index, neighbors[0].index),
            structure.get_distance(site_index, neighbors[1].index),
        ]
    )
    expected_bond_lengths = np.array(
        [
            float(1.36 + COVALENT_RADII[str(neighbors[0].site.specie)] - 0.76),
            float(1.36 + COVALENT_RADII[str(neighbors[1].site.specie)] - 0.76),
        ]
    )
    if str(neighbors[0].site.specie) == "H":
        num_h += 1
    if str(neighbors[1].site.specie) == "H":
        num_h += 1

    if (np.abs(180 - angle) < tolerance) or (np.abs(0 - angle) < tolerance):
        return False

    image_shift = np.zeros(3)
    for nn in connected_sites_a:
        if nn.index == site_index:
            image_shift = structure[site_index].coords - nn.site.coords
    dihedral_a = 0
    for nn in connected_sites_a:
        if nn.index != site_index:
            dihedral = vector_angle_to_plane(
                neighbors[1].site.coords,
                structure[site_index].coords,
                neighbors[0].site.coords,
                image_shift + nn.site.coords,
            )
            if dihedral > dihedral_a:
                dihedral_a = dihedral
    for nn in connected_sites_a:
        if nn.index == site_index:
            image_shift = structure[site_index].coords - nn.site.coords
    dihedral_b = 0
    for nn in connected_sites_b:
        if nn.index != site_index:
            dihedral = vector_angle_to_plane(
                neighbors[0].site.coords,
                structure[site_index].coords,
                neighbors[1].site.coords,
                image_shift + nn.site.coords,
            )
            if dihedral > dihedral_b:
                dihedral_b = dihedral
    mean_dihedral = np.min(np.abs([dihedral_a, dihedral_b]))
    if (np.abs(mean_dihedral - 180) < tolerance) or (np.abs(mean_dihedral - 0) < tolerance):
        if num_h == 1 and (bond_lengths[0] < expected_bond_lengths[0] and bond_lengths[1] < expected_bond_lengths[1]):
            return False
        if num_h == 0 and (bond_lengths[0] < expected_bond_lengths[0] or bond_lengths[1] < expected_bond_lengths[1]):
            return False
        if num_h == 0 and (bond_lengths[0] < (expected_bond_lengths[0] + 0.1) and bond_lengths[1] < (expected_bond_lengths[1] + 0.1)):
            return False
        return True
    return True


# ---------------------------------------------------------------------------
# Fused ring
# ---------------------------------------------------------------------------
def fused_ring_indices(structure, graph, cycles=None) -> list[int]:
    """Port of Fusedring_Check._get_fused_ring."""
    n_sum: list[int] = []
    cycles = _clean_cycles(graph) if cycles is None else cycles
    for site_index in _indices(structure, "N"):
        cn = _cn(graph, site_index)
        neighbors = _connected(graph, site_index)
        cm = _metal_neighbor_count(graph, site_index)
        if cn - cm != 2:
            continue
        non_metals = _non_metal_neighbors(neighbors)
        neighbors_index = [site_index] + [n.index for n in non_metals]
        n_5_ring = False
        for cycle in cycles:
            if len(cycle) == 5 and site_index in cycle and all(not _is_metal(structure[r]) for r in cycle):
                n_5_ring = True
                break
        if n_5_ring:
            for cycle in cycles:
                if (5 < len(cycle) < 10) and all(i in cycle for i in neighbors_index) and all(
                    not _is_metal(structure[r]) for r in cycle
                ):
                    n_sum.append(site_index)
    return n_sum


# ---------------------------------------------------------------------------
# Positive charge from linkers
# ---------------------------------------------------------------------------
def positive_charge_indices(structure, graph, cycles=None) -> list[int]:
    """Port of Positive_charge_Check._get_overcoordinated_nitrogen."""
    flagged: list[int] = []
    n_jump: list[int] = []
    cycles = _clean_cycles(graph) if cycles is None else cycles
    cycles16_by_node = _nonmetal_cycles_by_node(structure, cycles, 16)
    for site_index in _indices(structure, "N"):
        if site_index in n_jump:
            continue
        cn = _cn(graph, site_index)
        neighbors = _connected(graph, site_index)
        cm = _metal_neighbor_count(graph, site_index)
        if cn - cm == 4:
            flagged.append(site_index)
            continue
        for cycle in cycles16_by_node[site_index]:
            n_possible_jump: list[int] = []
            num_n = num_c = num_n3 = 0
            for ring in cycle:
                if _species(structure, ring) == "N":
                    n_possible_jump.append(ring)
                    num_n += 1
                    n_cn = _cn(graph, ring)
                    n_cm = _metal_neighbor_count(graph, ring)
                    if n_cn - n_cm == 3:
                        num_n3 += 1
                if _species(structure, ring) == "C":
                    num_c += 1
            if (num_n == 4) and (num_n3 == 4) and (num_c == 12):
                flagged.append(site_index)
                flagged.append(site_index)
            if (num_n == 4) and (num_n3 == 3) and (num_c == 12):
                flagged.append(site_index)
            if (num_n == 4) and (num_c == 12):
                n_jump.extend(n_possible_jump)
                break

    for site_index in _indices(structure, "O"):
        cn = _cn(graph, site_index)
        cm = _metal_neighbor_count(graph, site_index)
        if cn - cm == 3:
            flagged.append(site_index)
    for site_index in _indices(structure, "Sb"):
        flagged.extend([site_index] * 3)
    for site_index in _indices(structure, "Ge"):
        flagged.extend([site_index] * 4)
    return flagged


# ---------------------------------------------------------------------------
# Negative charge from linkers
# ---------------------------------------------------------------------------
def _neg_halogen(structure, graph) -> list[int]:
    halogen_sum: list[int] = []
    f_jump: list[int] = []
    halogen_indices = [i for i, s in enumerate(structure) if _is_halogen(s)]
    for site_index in halogen_indices:
        if site_index in f_jump:
            continue
        cn = _cn(graph, site_index)
        neighbors = _connected(graph, site_index)
        cm = _metal_neighbor_count(graph, site_index)
        non_metal = _non_metal_neighbor(neighbors)
        if cn - cm == 0:
            halogen_sum.append(site_index)
            continue
        if non_metal is None:
            continue
        cf = _num_neighbor_halogen(_connected(graph, non_metal.index))
        specie = str(non_metal.site.specie)
        if specie == "B" and cf == 4:
            halogen_sum.append(site_index)
            f_jump.extend(n.index for n in _connected(graph, non_metal.index))
        elif specie == "Si" and cf == 5:
            halogen_sum.append(site_index)
            f_jump.extend(n.index for n in _connected(graph, non_metal.index))
        elif specie == "Si" and cf == 6:
            halogen_sum.extend([site_index, site_index])
            f_jump.extend(n.index for n in _connected(graph, non_metal.index))
        elif specie == "P" and cf == 6:
            halogen_sum.append(site_index)
            f_jump.extend(n.index for n in _connected(graph, non_metal.index))
    return halogen_sum


def _neg_oxygen(structure, graph) -> list[int]:
    o_sum: list[int] = []
    o_jump: list[int] = []
    for site_index in _indices(structure, "O"):
        if site_index in o_jump:
            continue
        cn = _cn(graph, site_index)
        neighbors = _connected(graph, site_index)
        cm = _metal_neighbor_count(graph, site_index)
        if cn - cm == 0:
            o_sum.extend([site_index, site_index])
            continue
        if cn - cm == 2:
            continue
        if cn - cm != 1:
            continue
        non_metal = _non_metal_neighbor(neighbors)
        specie = str(non_metal.site.specie)
        if specie == "H":
            o_sum.append(site_index)
        elif specie == "C":
            neighbor_c = non_metal.index
            dis_c_o = structure.get_distance(neighbor_c, site_index)
            num_o = n_o = num_n = 0
            o_possible_jump = None
            for cnb in _connected(graph, neighbor_c):
                ncn = _cn(graph, cnb.index)
                ncm = _metal_neighbor_count(graph, cnb.index)
                sp = str(cnb.site.specie)
                if sp == "N":
                    num_n += 1
                if sp == "O":
                    n_o += 1
                if sp == "O" and ncn - ncm == 1:
                    num_o += 1
                    if cnb.index != site_index:
                        o_possible_jump = cnb.index
            if num_o == 2:
                o_sum.append(site_index)
                if o_possible_jump is not None:
                    o_jump.append(o_possible_jump)
            if n_o == 1:
                if dis_c_o > 1.315:
                    o_sum.append(site_index)
                elif dis_c_o <= 1.315 and num_n == 1:
                    o_sum.append(site_index)
        elif specie in ("N", "S", "P", "Cl", "Br", "I"):
            center = non_metal.index
            center_cn = _cn(graph, center)
            center_metal = _metal_neighbor_count(graph, center)
            o_attached = [n for n in _connected(graph, center) if str(n.site.specie) == "O"]
            num_o = len(o_attached)
            base = {"N": 1, "S": 2, "P": 3, "Cl": 1, "Br": 1, "I": 1}[specie]
            negative_o = base - (center_cn - center_metal - num_o)
            if negative_o > 0:
                for o_atom in o_attached:
                    o_cn = _cn(graph, o_atom.index)
                    o_cm = _metal_neighbor_count(graph, o_atom.index)
                    if o_atom.index != site_index:
                        o_jump.append(o_atom.index)
                    if o_cn - o_cm != 1:
                        negative_o -= 1
                # Charge contribution depends on the central element.
                if specie in ("N", "Cl", "Br", "I"):
                    if negative_o == 1:
                        o_sum.append(site_index)
                elif specie == "S":
                    if negative_o == 2:
                        o_sum.extend([site_index, site_index])
                    elif negative_o == 1:
                        o_sum.append(site_index)
                elif specie == "P":
                    if negative_o == 3:
                        o_sum.extend([site_index, site_index, site_index])
                    elif negative_o == 2:
                        o_sum.extend([site_index, site_index])
                    elif negative_o == 1:
                        o_sum.append(site_index)
        else:
            o_sum.append(site_index)
    return o_sum


def _neg_sulfur(structure, graph) -> list[int]:
    s_sum: list[int] = []
    s_jump: list[int] = []
    for site_index in _indices(structure, "S"):
        if site_index in s_jump:
            continue
        cn = _cn(graph, site_index)
        neighbors = _connected(graph, site_index)
        cm = _metal_neighbor_count(graph, site_index)
        if cn - cm == 0:
            s_sum.extend([site_index, site_index])
            continue
        if cn - cm == 2:
            continue
        if cn - cm != 1:
            continue
        non_metal = _non_metal_neighbor(neighbors)
        specie = str(non_metal.site.specie)
        if specie == "H":
            s_sum.append(site_index)
        elif specie == "C":
            neighbor_c = non_metal.index
            dis_c_s = structure.get_distance(neighbor_c, site_index)
            num_s = 1
            negative_s = 0
            num_n = 0
            s_possible_jump = None
            for cnb in _connected(graph, neighbor_c):
                ncn = _cn(graph, cnb.index)
                ncm = _metal_neighbor_count(graph, cnb.index)
                sp = str(cnb.site.specie)
                if sp == "S" and cnb.index != site_index:
                    num_s += 1
                    s_possible_jump = cnb.index
                if sp == "N":
                    num_n += 1
                if ncn - ncm == 1:
                    negative_s += 1
            if num_s == 1:
                if dis_c_s > 1.66:
                    s_sum.append(site_index)
                elif dis_c_s <= 1.66 and num_n == 1:
                    s_sum.append(site_index)
            if num_s == 2 and negative_s >= 1:
                s_sum.append(site_index)
                if s_possible_jump is not None:
                    s_jump.append(s_possible_jump)
        else:
            s_sum.append(site_index)
    return s_sum


def _neg_nitrogen(structure, graph, cycles=None) -> list[int]:
    n_sum: list[int] = []
    n_jump: list[int] = []
    cycles = _clean_cycles(graph) if cycles is None else cycles
    for site_index in _indices(structure, "N"):
        if site_index in n_jump:
            continue
        cn = _cn(graph, site_index)
        neighbors = _connected(graph, site_index)
        cm = _metal_neighbor_count(graph, site_index)
        if cn - cm == 1:
            non_metal = _non_metal_neighbor(neighbors)
            specie = str(non_metal.site.specie)
            if specie == "C":
                neighbor_c = non_metal.index
                cn_c = _cn(graph, neighbor_c)
                cm_c = _metal_neighbor_count(graph, neighbor_c)
                if cn_c - cm_c == 1:
                    n_sum.append(site_index)
            if specie == "N":
                neighbor_n = non_metal.index
                nn_neighbors = _connected(graph, neighbor_n)
                cn_n = _cn(graph, neighbor_n)
                if cn_n == 2 and str(nn_neighbors[0].site.specie) == "N" and str(nn_neighbors[1].site.specie) == "N":
                    angle = get_angle_between_site_and_neighbors(non_metal.site, nn_neighbors)
                    if np.abs(180 - angle) < 25 or np.abs(0 - angle) < 25:
                        n_sum.append(site_index)
                        n_jump.append(nn_neighbors[0].index)
                        n_jump.append(nn_neighbors[1].index)
        if cn - cm == 2:
            non_metals = _non_metal_neighbors(neighbors)
            undercoordinated = _guess_underbound_nitrogen_cn2(
                structure,
                site_index,
                non_metals,
                _connected(graph, non_metals[0].index),
                _connected(graph, non_metals[1].index),
                30,
            )
            if undercoordinated:
                n_sum.append(site_index)
                continue
            for cycle in cycles:
                n_possible_jump: list[int] = []
                if len(cycle) == 16 and site_index in cycle and all(not _is_metal(structure[r]) for r in cycle):
                    num_n = num_c = num_n3 = 0
                    for ring in cycle:
                        if _species(structure, ring) == "N":
                            n_possible_jump.append(ring)
                            num_n += 1
                            r_cn = _cn(graph, ring)
                            r_cm = _metal_neighbor_count(graph, ring)
                            if r_cn - r_cm == 3:
                                num_n3 += 1
                        if _species(structure, ring) == "C":
                            num_c += 1
                    if (num_n == 4) and (num_n3 == 0) and (num_c == 12):
                        n_sum.extend([site_index, site_index])
                    if (num_n == 4) and (num_n3 == 1) and (num_c == 12):
                        n_sum.append(site_index)
                    if (num_n == 4) and (num_c == 12):
                        n_jump.extend(n_possible_jump)
                        break
            if site_index not in n_jump:
                for cycle in cycles:
                    cn3_n = 0
                    if len(cycle) == 5 and site_index in cycle and all(not _is_metal(structure[r]) for r in cycle):
                        for ring in cycle:
                            ring_cn = len(_non_metal_neighbors(_connected(graph, ring)))
                            if _species(structure, ring) == "N":
                                n_jump.append(ring)
                            if _species(structure, ring) == "N" and ring_cn == 3:
                                cn3_n += 1
                            if _species(structure, ring) in ("O", "S"):
                                cn3_n += 1
                        if cn3_n in (0, 2, 4):
                            n_sum.append(site_index)
    return n_sum


def negative_charge_indices(structure, graph, cycles=None) -> list[int]:
    """Port of Negative_charge_Check._run_check (halogen + O + S + N)."""
    return (
        _neg_halogen(structure, graph)
        + _neg_oxygen(structure, graph)
        + _neg_sulfur(structure, graph)
        + _neg_nitrogen(structure, graph, cycles=cycles)
    )


# ---------------------------------------------------------------------------
# Open metal sites
# ---------------------------------------------------------------------------
def _check_if_open(lsop, is_open, weights, threshold: float = 0.5):
    if lsop is None:
        return None
    if is_open is None:
        return False
    lsop = np.array(lsop) * np.array(weights)
    open_contributions = lsop[is_open].sum()
    close_contributions = lsop.sum() - open_contributions
    return open_contributions / (open_contributions + close_contributions) > threshold


@lru_cache(maxsize=None)
def _lsop(names: tuple[str, ...]):
    from pymatgen.analysis.local_env import LocalStructOrderParams

    return LocalStructOrderParams(list(names))


def _ops_for_site(structure, graph, site_index):
    cn = _cn(graph, site_index)
    if cn not in OP_DEF:
        # Mirror the reference: CN<=3 -> open, CN>8 -> undefined (None).
        return cn, None, None, None, None
    names = tuple(OP_DEF[cn]["names"])
    is_open = OP_DEF[cn]["open"]
    weights = OP_DEF[cn]["weights"]
    return cn, names, _lsop(names).get_order_parameters(structure, site_index), is_open, weights


def _is_open_metal_site(structure, graph, site_index) -> bool:
    cn = _cn(graph, site_index)
    if cn <= 3:
        return True
    if cn > 8:
        return False
    _, _, lsop, is_open, weights = _ops_for_site(structure, graph, site_index)
    return bool(_check_if_open(lsop, is_open, weights))


def oms_indices(structure, graph) -> list[int]:
    """Port of MOFOMS.check_oms."""
    return [i for i, site in enumerate(structure) if _is_metal(site) and _is_open_metal_site(structure, graph, i)]


def has_oms(structure, graph) -> bool:
    return any(_is_open_metal_site(structure, graph, i) for i, site in enumerate(structure) if _is_metal(site))


# ---------------------------------------------------------------------------
# Structure-level entry points (build the graph the reference way).
# ---------------------------------------------------------------------------
def possible_charged_fused_ring_from_structure(structure, method: str = "vesta") -> bool:
    graph = build_structure_graph(structure, method)
    return len(fused_ring_indices(structure, graph)) > 0


def positive_charge_from_linkers_from_structure(structure, method: str = "vesta") -> int:
    graph = build_structure_graph(structure, method)
    return len(positive_charge_indices(structure, graph))


def negative_charge_from_linkers_from_structure(structure, method: str = "vesta") -> int:
    graph = build_structure_graph(structure, method)
    return len(negative_charge_indices(structure, graph))


def has_oms_from_structure(structure, method: str = "vesta") -> bool:
    graph = build_structure_graph(structure, method)
    return has_oms(structure, graph)


def oms_indices_from_structure(structure, method: str = "vesta") -> list[int]:
    graph = build_structure_graph(structure, method)
    return oms_indices(structure, graph)
