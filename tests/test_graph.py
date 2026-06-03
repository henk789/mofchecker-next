import random
from types import SimpleNamespace

import pytest

from mofchecker_next._rust import connected_components
from mofchecker_next.checks.graph import (
    check_nonperiodic_components,
    check_overcoordinated_by_degree,
    connected_components_from_edges,
    connected_components_py,
    false_oxo_indices_by_graph,
    geometrically_exposed_metal_indices_from_graph,
    get_open_angle_from_coords,
    is_3d_connected_graph_from_structure,
    nonperiodic_component_indices,
    node_degrees_from_edges,
    overcoordinated_indices_by_degree,
)


def test_connected_components_empty_graph():
    assert connected_components(0, []) == []


def test_connected_components_isolated_atoms():
    assert connected_components(3, []) == [[0], [1], [2]]


def test_connected_components_multiple_components_stable_ordering():
    edges = [(4, 5), (1, 2), (2, 3)]
    assert connected_components(6, edges) == [[0], [1, 2, 3], [4, 5]]


def test_connected_components_ignores_self_and_duplicate_edges():
    edges = [(0, 0), (0, 1), (1, 0), (1, 2)]
    assert connected_components(3, edges) == [[0, 1, 2]]


def test_connected_components_matches_python_randomized():
    rng = random.Random(314159)
    for _ in range(100):
        n_atoms = rng.randint(0, 20)
        edges = []
        for _ in range(rng.randint(0, 40)):
            if n_atoms == 0:
                break
            edges.append((rng.randrange(n_atoms), rng.randrange(n_atoms)))
        assert connected_components(n_atoms, edges) == connected_components_py(n_atoms, edges)


def test_connected_components_rejects_out_of_bounds_edges():
    with pytest.raises(Exception):
        connected_components(2, [(0, 2)])


def test_connected_components_from_edges_uses_rust_kernel():
    assert connected_components_from_edges(4, [(0, 1), (2, 3)]) == [[0, 1], [2, 3]]


def test_nonperiodic_component_indices_use_edge_images():
    edges = [(0, 1), (2, 3), (3, 2)]
    edge_images = [(1, 0, 0), (0, 0, 0), (0, 0, 0)]
    assert nonperiodic_component_indices(4, edges, edge_images) == [[2, 3]]


def test_nonperiodic_component_indices_include_isolated_atoms():
    assert nonperiodic_component_indices(3, [(0, 1)], [(1, 0, 0)]) == [[2]]


def test_nonperiodic_component_indices_validate_image_count():
    with pytest.raises(ValueError):
        nonperiodic_component_indices(2, [(0, 1)], [])


def test_check_nonperiodic_components_returns_diagnostics():
    diagnostics = check_nonperiodic_components(3, [(0, 1)], [(1, 0, 0)])
    assert len(diagnostics) == 1
    assert diagnostics[0].check == "floating_component"
    assert diagnostics[0].atoms[0].index == 2
    assert diagnostics[0].values == {"component_size": 1}


def test_nonperiodic_component_multiple_components():
    edges = [(0, 1), (2, 3), (4, 5)]
    edge_images = [(1, 0, 0), (0, 0, 0), (0, 0, 0)]
    assert nonperiodic_component_indices(6, edges, edge_images) == [[2, 3], [4, 5]]


def test_node_degrees_from_edges_counts_unique_neighbors():
    assert node_degrees_from_edges(4, [(0, 1), (1, 0), (1, 2), (2, 2)]) == [1, 2, 1, 0]


def test_overcoordinated_indices_by_degree_flags_target_only():
    atomic_numbers = [6, 1, 1, 1, 1, 8]
    edges = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)]
    assert overcoordinated_indices_by_degree(atomic_numbers, edges, 6, 4) == [0]
    assert overcoordinated_indices_by_degree(atomic_numbers, edges, 8, 2) == []


def test_overcoordinated_nitrogen_positive_and_metal_exclusion():
    atomic_numbers = [7, 1, 1, 1, 1, 1, 26]
    edges = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)]
    assert overcoordinated_indices_by_degree(atomic_numbers, edges, 7, 4) == [0]
    edges_with_metal = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 6)]
    assert overcoordinated_indices_by_degree(atomic_numbers, edges_with_metal, 7, 4, [26]) == []


def test_overcoordinated_indices_by_degree_applies_neighbor_exclusions():
    atomic_numbers = [6, 1, 1, 1, 1, 5]
    edges = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)]
    assert overcoordinated_indices_by_degree(atomic_numbers, edges, 6, 4) == [0]
    assert overcoordinated_indices_by_degree(atomic_numbers, edges, 6, 4, [5]) == []


def test_false_oxo_indices_by_graph_flags_disallowed_terminal_oxo_metal():
    atomic_symbols = ["Fe", "O", "C"]
    edges = [(0, 1)]
    assert false_oxo_indices_by_graph(atomic_symbols, edges, {"Fe"}) == [0]


def test_false_oxo_indices_by_graph_ignores_allowed_or_nonterminal_oxo():
    atomic_symbols = ["V", "O", "Fe", "O", "C"]
    edges = [(0, 1), (2, 3), (3, 4)]
    assert false_oxo_indices_by_graph(atomic_symbols, edges, {"V", "Fe"}) == []


def test_get_open_angle_from_coords_returns_reference_coplanar_fallback():
    coords = [[0, 0, 0], [1, 0, 0], [-1, 0, 0]]
    assert get_open_angle_from_coords(coords, ["Fe", "O", "O"]) == 180


def test_geometrically_exposed_metal_indices_from_graph_flags_open_low_cn_metal():
    class FakeStructure:
        def __init__(self, sites, cart_coords):
            self._sites = sites
            self.cart_coords = cart_coords

        def __iter__(self):
            return iter(self._sites)

        def __getitem__(self, index):
            return self._sites[index]

    sites = [SimpleNamespace(specie="Fe"), SimpleNamespace(specie="O"), SimpleNamespace(specie="O")]
    structure = FakeStructure(sites, [[0, 0, 0], [1, 0, 0], [-1, 0, 0]])
    graph = SimpleNamespace(structure=structure)
    neighbors = {
        0: [SimpleNamespace(index=1, site=SimpleNamespace(specie="O", coords=[1, 0, 0])), SimpleNamespace(index=2, site=SimpleNamespace(specie="O", coords=[-1, 0, 0]))],
        1: [SimpleNamespace(index=0, site=SimpleNamespace(specie="Fe", coords=[0, 0, 0]))],
        2: [SimpleNamespace(index=0, site=SimpleNamespace(specie="Fe", coords=[0, 0, 0]))],
    }
    graph.get_connected_sites = lambda index: neighbors[index]

    assert geometrically_exposed_metal_indices_from_graph(graph, {"Fe"}) == [0]


def test_check_overcoordinated_by_degree_returns_diagnostics():
    atomic_numbers = [7, 1, 1, 1, 1, 1]
    edges = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)]
    diagnostics = check_overcoordinated_by_degree(
        atomic_numbers,
        edges,
        target_atomic_number=7,
        max_degree=4,
        check_name="overcoordinated_nitrogen",
    )
    assert len(diagnostics) == 1
    assert diagnostics[0].atoms[0].index == 0
    assert diagnostics[0].values == {"coordination_number": 5, "max_allowed": 4}


def test_is_3d_connected_graph_from_structure_returns_bool():
    from pymatgen.core import Lattice, Structure

    structure = Structure(Lattice.cubic(10.0), ["Br"], [[0, 0, 0]])
    assert is_3d_connected_graph_from_structure(structure) is False
