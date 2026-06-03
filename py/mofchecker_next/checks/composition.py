from __future__ import annotations


def species_symbols(structure) -> list[str]:
    return [str(species) for species in structure.species]


def element_indices(structure, symbol: str) -> list[int]:
    return [index for index, species in enumerate(species_symbols(structure)) if species == symbol]


def metal_indices(structure, metal_symbols) -> list[int]:
    metals = {str(symbol) for symbol in metal_symbols}
    return [index for index, site in enumerate(structure) if str(site.specie) in metals]


def has_element(structure, symbol: str) -> bool:
    return len(element_indices(structure, symbol)) > 0


def has_metal(structure, metal_symbols) -> bool:
    return len(metal_indices(structure, metal_symbols)) > 0


def metal_number(structure, metal_symbols) -> int:
    """Match MOFChecker 2.0's metal number convention: metals minus Sb."""
    return len(metal_indices(structure, metal_symbols)) - len(element_indices(structure, "Sb"))


def simple_global_diagnostics(structure, metal_symbols) -> dict:
    return {
        "has_carbon": has_element(structure, "C"),
        "has_hydrogen": has_element(structure, "H"),
        "has_nitrogen": has_element(structure, "N"),
        "has_metal": has_metal(structure, metal_symbols),
        "metal_number": metal_number(structure, metal_symbols),
    }
