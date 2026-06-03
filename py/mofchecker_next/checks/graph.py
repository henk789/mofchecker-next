from __future__ import annotations

import math

import numpy as np
from element_coder.encode import encode_many
from libconeangle import cone_angle

from mofchecker_next.diagnostics import AtomRef, Diagnostic


NO_TERMINAL_OXO_SYMBOLS = {
    "Li",
    "Na",
    "K",
    "Rb",
    "Cs",
    "Fr",
    "Be",
    "Mg",
    "Ca",
    "Sr",
    "Ba",
    "Ra",
    "Sc",
    "Y",
    "La",
    "Ac",
    "Ti",
    "Zr",
    "Hf",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Ag",
    "Zn",
    "Cd",
    "Al",
    "Ga",
    "In",
    "Tl",
}


def build_structure_graph(structure, method: str = "vesta"):
    """Build a StructureGraph using the MOFChecker 2.0 graph helper."""
    from structuregraph_helpers.create import get_structure_graph

    return get_structure_graph(structure, method)


def structure_graph_edges_and_images(structure_graph):
    """Extract explicit edge endpoints and `to_jimage` vectors from a StructureGraph."""
    edges = []
    edge_images = []
    for u, v, data in structure_graph.graph.edges(data=True):
        edges.append((int(u), int(v)))
        edge_images.append(tuple(int(x) for x in data.get("to_jimage", (0, 0, 0))))
    return edges, edge_images


def structure_graph_components_as_molecules(structure_graph):
    """Call the same helper MOFChecker 2.0 uses for floating-solvent indices."""
    from structuregraph_helpers.subgraph import get_subgraphs_as_molecules

    _, _, indices, _, _ = get_subgraphs_as_molecules(structure_graph, return_unique=False)
    return indices


def floating_solvent_indices_from_structure(structure, method: str = "vesta"):
    """Return MOFChecker-compatible floating-solvent index lists."""
    return structure_graph_components_as_molecules(build_structure_graph(structure, method))


def is_3d_connected_graph_from_structure(structure, method: str = "vesta") -> bool:
    """Match MOFChecker 2.0's 3D connected graph descriptor."""
    from pymatgen.analysis.dimensionality import get_dimensionality_larsen

    return bool(get_dimensionality_larsen(build_structure_graph(structure, method)) == 3)


def connected_components_py(n_atoms: int, edges):
    """Pure-Python reference connected components over explicit atom-index edges."""
    adjacency = [[] for _ in range(n_atoms)]
    for a, b in edges:
        a = int(a)
        b = int(b)
        if not (0 <= a < n_atoms and 0 <= b < n_atoms):
            raise IndexError("edge endpoint out of bounds")
        if a == b:
            continue
        adjacency[a].append(b)
        adjacency[b].append(a)
    adjacency = [sorted(set(neighbors)) for neighbors in adjacency]

    seen = [False] * n_atoms
    components = []
    for start in range(n_atoms):
        if seen[start]:
            continue
        stack = [start]
        seen[start] = True
        component = []
        while stack:
            node = stack.pop()
            component.append(node)
            for neighbor in reversed(adjacency[node]):
                if not seen[neighbor]:
                    seen[neighbor] = True
                    stack.append(neighbor)
        components.append(sorted(component))
    return sorted(components, key=lambda component: component[0])


def connected_components_from_edges(n_atoms: int, edges):
    """Run the Rust connected-components kernel on explicit Python-built edges."""
    from mofchecker_next._rust import connected_components

    return connected_components(int(n_atoms), [(int(a), int(b)) for a, b in edges])


def node_degrees_from_edges(n_atoms: int, edges) -> list[int]:
    """Return unique graph degree for each node using explicit edges."""
    from mofchecker_next._rust import node_degrees

    return node_degrees(int(n_atoms), [(int(a), int(b)) for a, b in edges])


def overcoordinated_indices_by_degree(
    atomic_numbers,
    edges,
    target_atomic_number: int,
    max_degree: int,
    excluded_neighbor_atomic_numbers=(),
) -> list[int]:
    """Flag target atoms with graph degree above a threshold.

    Neighbor exclusions are caller-supplied so Rust does not infer metals,
    boron chemistry, or bond perception rules.
    """
    atomic_numbers = [int(number) for number in atomic_numbers]
    edges = [(int(a), int(b)) for a, b in edges]
    excluded = {int(number) for number in excluded_neighbor_atomic_numbers}
    degrees = node_degrees_from_edges(len(atomic_numbers), edges)
    neighbors = [set() for _ in atomic_numbers]
    for a, b in edges:
        if a == b:
            continue
        neighbors[a].add(b)
        neighbors[b].add(a)

    flagged = []
    for index, atomic_number in enumerate(atomic_numbers):
        if atomic_number != int(target_atomic_number):
            continue
        if degrees[index] <= int(max_degree):
            continue
        if any(atomic_numbers[neighbor] in excluded for neighbor in neighbors[index]):
            continue
        flagged.append(index)
    return flagged


def false_oxo_indices_by_graph(
    atomic_symbols,
    edges,
    metal_symbols,
    no_terminal_oxo_symbols=NO_TERMINAL_OXO_SYMBOLS,
) -> list[int]:
    """Return metals with a disallowed terminal O neighbor on an explicit graph."""
    atomic_symbols = [str(symbol) for symbol in atomic_symbols]
    metal_symbols = {str(symbol) for symbol in metal_symbols}
    no_terminal_oxo_symbols = {str(symbol) for symbol in no_terminal_oxo_symbols}
    neighbors = [set() for _ in atomic_symbols]
    for a, b in edges:
        a = int(a)
        b = int(b)
        if a == b:
            continue
        neighbors[a].add(b)
        neighbors[b].add(a)

    flagged = []
    for site_index, symbol in enumerate(atomic_symbols):
        if symbol not in metal_symbols or symbol not in no_terminal_oxo_symbols:
            continue
        for neighbor_index in neighbors[site_index]:
            if atomic_symbols[neighbor_index] == "O" and len(neighbors[neighbor_index]) == 1:
                flagged.append(site_index)
    return flagged


def get_open_angle_from_coords(coords, species) -> float:
    """Return MOFChecker 2.0's open angle, defined as 360 minus cone angle."""
    coords = np.array(coords, dtype=np.float64)
    species = [str(symbol) for symbol in species]
    encodings = np.array(encode_many(species, "van_der_waals_radius"), dtype=np.float64)
    try:
        angle, _, _ = cone_angle(coords, encodings, 0)
        return 360 - angle
    except ValueError:
        centered = np.unique(np.array(coords), axis=0)
        centered -= centered.mean(axis=0)
        if np.linalg.matrix_rank(centered) <= 2:
            return 180
        return np.nan


def get_open_angle_from_graph(structure_graph, site_index: int) -> float:
    coords = [structure_graph.structure.cart_coords[int(site_index)]]
    species = [str(structure_graph.structure[int(site_index)].specie)]
    for neighbor in structure_graph.get_connected_sites(int(site_index)):
        coords.append(neighbor.site.coords)
        species.append(str(neighbor.site.specie))
    return get_open_angle_from_coords(coords, species)


def geometrically_exposed_metal_indices_from_graph(
    structure_graph,
    metal_symbols,
    threshold: float = 150.0,
) -> list[int]:
    """Return metal indices with open angle above threshold and graph CN below 6."""
    metal_symbols = {str(symbol) for symbol in metal_symbols}
    flagged = []
    for site_index, site in enumerate(structure_graph.structure):
        if str(site.specie) not in metal_symbols:
            continue
        angle = get_open_angle_from_graph(structure_graph, site_index)
        if angle > threshold and len(structure_graph.get_connected_sites(site_index)) < 6:
            flagged.append(site_index)
    return flagged


def check_overcoordinated_by_degree(
    atomic_numbers,
    edges,
    target_atomic_number: int,
    max_degree: int,
    check_name: str,
    excluded_neighbor_atomic_numbers=(),
) -> list[Diagnostic]:
    """Return diagnostics for a simple explicit-graph overcoordination rule."""
    degrees = node_degrees_from_edges(len(atomic_numbers), edges)
    diagnostics = []
    for index in overcoordinated_indices_by_degree(
        atomic_numbers,
        edges,
        target_atomic_number,
        max_degree,
        excluded_neighbor_atomic_numbers,
    ):
        diagnostics.append(
            Diagnostic(
                check=check_name,
                severity="error",
                message="Atom has more graph neighbors than allowed",
                atoms=[AtomRef(index=index)],
                values={"coordination_number": degrees[index], "max_allowed": int(max_degree)},
            )
        )
    return diagnostics


def nonperiodic_component_indices(n_atoms: int, edges, edge_images=None):
    """Return components whose internal edges do not cross periodic boundaries.

    This is a small explicit-edge helper for diagnostics. It is not a full port
    of MOFChecker 2.0's floating-solvent algorithm, which operates on a 3x3x3
    supercell through ``structuregraph_helpers``.
    """
    edges = [(int(a), int(b)) for a, b in edges]
    if edge_images is None:
        edge_images = [(0, 0, 0)] * len(edges)
    edge_images = [tuple(int(v) for v in image) for image in edge_images]
    if len(edges) != len(edge_images):
        raise ValueError("edges and edge_images must have the same length")

    components = connected_components_from_edges(n_atoms, edges)
    node_to_component = {}
    for component_index, component in enumerate(components):
        for node in component:
            node_to_component[node] = component_index

    crosses_boundary = [False] * len(components)
    for (a, b), image in zip(edges, edge_images):
        if a == b:
            continue
        component_index = node_to_component[a]
        if component_index == node_to_component[b] and image != (0, 0, 0):
            crosses_boundary[component_index] = True

    return [component for component, crosses in zip(components, crosses_boundary) if not crosses]


def check_nonperiodic_components(n_atoms: int, edges, edge_images=None) -> list[Diagnostic]:
    """Return diagnostics for finite/non-periodic connected components."""
    diagnostics = []
    for component in nonperiodic_component_indices(n_atoms, edges, edge_images):
        diagnostics.append(
            Diagnostic(
                check="floating_component",
                severity="warning",
                message="Non-periodic connected component",
                atoms=[AtomRef(index=index) for index in component],
                values={"component_size": len(component)},
            )
        )
    return diagnostics


def overcoordinated_carbon_indices_from_structure(structure, metal_symbols, method: str = "vesta"):
    graph = build_structure_graph(structure, method)
    edges, _ = structure_graph_edges_and_images(graph)
    atomic_numbers = [int(site.specie.Z) for site in graph.structure]
    excluded = {5}
    excluded.update(int(site.specie.Z) for site in graph.structure if str(site.specie) in set(metal_symbols))
    return overcoordinated_indices_by_degree(atomic_numbers, edges, 6, 4, excluded)


def overcoordinated_nitrogen_indices_from_structure(structure, metal_symbols, method: str = "vesta"):
    graph = build_structure_graph(structure, method)
    edges, _ = structure_graph_edges_and_images(graph)
    atomic_numbers = [int(site.specie.Z) for site in graph.structure]
    excluded = {int(site.specie.Z) for site in graph.structure if str(site.specie) in set(metal_symbols)}
    return overcoordinated_indices_by_degree(atomic_numbers, edges, 7, 4, excluded)


def false_oxo_indices_from_structure(
    structure,
    metal_symbols,
    method: str = "vesta",
    no_terminal_oxo_symbols=NO_TERMINAL_OXO_SYMBOLS,
) -> list[int]:
    """Return terminal oxo indices disallowed by MOFChecker 2.0's element list."""
    graph = build_structure_graph(structure, method)
    metal_symbols = {str(symbol) for symbol in metal_symbols}
    no_terminal_oxo_symbols = {str(symbol) for symbol in no_terminal_oxo_symbols}

    atomic_symbols = [str(site.specie) for site in structure]
    edges, _ = structure_graph_edges_and_images(graph)
    return false_oxo_indices_by_graph(atomic_symbols, edges, metal_symbols, no_terminal_oxo_symbols)


def geometrically_exposed_metal_indices_from_structure(
    structure,
    metal_symbols,
    method: str = "vesta",
) -> list[int]:
    graph = build_structure_graph(structure, method)
    return geometrically_exposed_metal_indices_from_graph(graph, metal_symbols)


def undercoordinated_carbon_indices_from_structure(
    structure,
    metal_symbols,
    covalent_radii_by_symbol,
    method: str = "vesta",
    tolerance: float = 165.0,
) -> list[int]:
    """Return undercoordinated carbon indices following MOFChecker 2.0 index logic.

    This intentionally matches only flagged indices. MOFChecker's candidate H
    positions are correction/healing scaffolding and are outside current scope.
    """
    graph = build_structure_graph(structure, method)
    c_indices = [index for index, site in enumerate(structure) if str(site.specie) == "C"]
    n_indices = {index for index, site in enumerate(structure) if str(site.specie) == "N"}
    metal_symbols = {str(symbol) for symbol in metal_symbols}
    radii = {str(symbol): float(radius) for symbol, radius in covalent_radii_by_symbol.items()}

    flagged = []
    for site_index in c_indices:
        neighbors = graph.get_connected_sites(site_index)
        cn = len(neighbors)
        if cn == 1:
            if neighbors[0].index not in n_indices:
                flagged.append(site_index)
            elif structure.get_distance(neighbors[0].index, site_index) > 1.2:
                flagged.append(site_index)

    for site_index in c_indices:
        neighbors = graph.get_connected_sites(site_index)
        cn = len(neighbors)
        if cn != 2:
            continue
        a = structure.get_distance(neighbors[0].index, site_index)
        b = structure.get_distance(neighbors[1].index, site_index)
        c = structure.get_distance(neighbors[0].index, neighbors[1].index)
        cos_angle = (a * a + b * b - c * c) / (2 * a * b)
        angle = math.degrees(math.acos(round(cos_angle, 6)))
        any_metal_neighbor = any(str(neighbor.site.specie) in metal_symbols for neighbor in neighbors)
        if any_metal_neighbor:
            if angle < tolerance - 15:
                flagged.append(site_index)
        elif angle < tolerance:
            flagged.append(site_index)

    return flagged


def undercoordinated_nitrogen_indices_from_structure(
    structure,
    metal_symbols,
    method: str = "vesta",
) -> list[int]:
    """Return undercoordinated nitrogen indices following MOFChecker 2.0 index logic.

    MOFChecker's CN=2 nitrogen branch currently only records candidate positions,
    not flagged nitrogen indices. This helper preserves that reference-visible
    behavior and returns indices from the CN-minus-metal-equals-1 branch.
    """
    graph = build_structure_graph(structure, method)
    metal_symbols = {str(symbol) for symbol in metal_symbols}
    n_indices = [index for index, site in enumerate(structure) if str(site.specie) == "N"]
    flagged = []
    for site_index in n_indices:
        neighbors = graph.get_connected_sites(site_index)
        metal_count = sum(1 for neighbor in neighbors if str(neighbor.site.specie) in metal_symbols)
        if len(neighbors) - metal_count != 1:
            continue
        non_metal_neighbors = [neighbor for neighbor in neighbors if str(neighbor.site.specie) not in metal_symbols]
        neighbor = non_metal_neighbors[0]
        neighbor_species = str(neighbor.site.specie)
        distance = structure.get_distance(neighbor.index, site_index)
        if neighbor_species in {"C", "N"} and distance < 1.25:
            continue
        if neighbor_species in {"C", "N"} and distance > 1.35:
            flagged.append(site_index)
        elif neighbor_species not in {"C", "N"} or not distance > 1.35:
            flagged.append(site_index)
    return flagged


def undercoordinated_rare_earth_indices_from_structure(structure, method: str = "vesta") -> list[int]:
    graph = build_structure_graph(structure, method)
    degrees = node_degrees_from_edges(len(structure), structure_graph_edges_and_images(graph)[0])
    return [index for index, site in enumerate(structure) if site.specie.is_rare_earth and degrees[index] < 4]


def undercoordinated_alkali_alkaline_indices_from_structure(structure, method: str = "vesta") -> list[int]:
    graph = build_structure_graph(structure, method)
    degrees = node_degrees_from_edges(len(structure), structure_graph_edges_and_images(graph)[0])
    return [
        index
        for index, site in enumerate(structure)
        if (site.specie.is_alkali or site.specie.is_alkaline) and degrees[index] < 4
    ]
