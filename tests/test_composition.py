from pymatgen.core import Lattice, Structure

from mofchecker_next.checks.composition import (
    element_indices,
    has_element,
    has_metal,
    metal_indices,
    metal_number,
    simple_global_diagnostics,
)


def test_element_indices_and_has_element():
    structure = Structure(Lattice.cubic(10.0), ["C", "H", "C"], [[0, 0, 0], [0.1, 0, 0], [0.2, 0, 0]])
    assert element_indices(structure, "C") == [0, 2]
    assert has_element(structure, "H") is True
    assert has_element(structure, "N") is False


def test_metal_indices_use_supplied_metal_set():
    structure = Structure(Lattice.cubic(10.0), ["Na", "C", "Sb"], [[0, 0, 0], [0.1, 0, 0], [0.2, 0, 0]])
    metal_symbols = {"Na", "Sb"}
    assert metal_indices(structure, metal_symbols) == [0, 2]
    assert has_metal(structure, metal_symbols) is True
    assert metal_number(structure, metal_symbols) == 1


def test_simple_global_diagnostics():
    structure = Structure(Lattice.cubic(10.0), ["Na", "C", "H"], [[0, 0, 0], [0.1, 0, 0], [0.2, 0, 0]])
    assert simple_global_diagnostics(structure, {"Na"}) == {
        "has_carbon": True,
        "has_hydrogen": True,
        "has_nitrogen": False,
        "has_metal": True,
        "metal_number": 1,
    }
