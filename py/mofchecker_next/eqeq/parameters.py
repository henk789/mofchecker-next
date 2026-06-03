"""Per-atom EQeq parameters from the vendored ionization/charge-center tables.

Faithful translation of the parameter setup in EQeq's ``main.cpp``
(``LoadIonizationData``, ``LoadChargeCenters`` and the ``X``/``J`` assignment in
``LoadCIFData``). The data files ``data/ionizationdata.dat`` and
``data/chargecenters.dat`` are vendored verbatim from EQeq (GPLv2 -- see the
LICENSE in this subpackage), so the derived electronegativity and hardness match
the reference exactly.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

_DATA_DIR = Path(__file__).resolve().parent / "data"
_IONIZATION_FILE = _DATA_DIR / "ionizationdata.dat"
_CHARGE_CENTER_FILE = _DATA_DIR / "chargecenters.dat"

# EQeq's two global parameters.
DEFAULT_LAMBDA = 1.2  # dielectric screening (eps_eff = 1.67)
DEFAULT_H_ELECTRON_AFFINITY = -2.0  # hI0; hydrogen's effective electron affinity

# Hydrogen constant from main.cpp (empirical 1st ionization energy).
_HI1 = 13.598

# EQeq's element ordering (enum StringAtomLabels, 0-based: H=0). The ionization
# table is read in this order and charge centers are looked up by it.
_ELEMENT_ORDER = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni "
    "Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I "
    "Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt "
    "Au Hg Tl Pb Bi Po"
).split()
_ELEMENT_INDEX = {symbol: i for i, symbol in enumerate(_ELEMENT_ORDER)}


@dataclass(frozen=True)
class AtomParameters:
    electronegativity: float  # X
    hardness: float  # J


def _parse_token(token: str) -> float | None:
    """Parse one ionization-data token; ``na``/``np`` mean unavailable."""
    token = token.strip()
    if token in ("na", "np", ""):
        return None
    if token.startswith("<0.5"):
        return 0.5
    return float(token)


@lru_cache(maxsize=1)
def _ionization_table() -> list[list[float]]:
    """Return ``ionizationPotential[Z] = [EA, IE1, ..., IE8]`` (9 entries)."""
    table: list[list[float]] = []
    for line in _IONIZATION_FILE.read_text().splitlines():
        fields = line.strip().split("\t")
        if len(fields) < 4:
            continue
        # fields: Z, label, status, EA, IE1..IE8
        values: list[float] = []
        for token in fields[3:12]:
            parsed = _parse_token(token)
            values.append(parsed if parsed is not None else 0.0)
        while len(values) < 9:
            values.append(0.0)
        table.append(values)
    return table


@lru_cache(maxsize=1)
def _charge_centers() -> dict[str, int]:
    centers: dict[str, int] = {}
    for line in _CHARGE_CENTER_FILE.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        centers[parts[0]] = int(parts[1])
    return centers


@lru_cache(maxsize=256)
def parameters_for(
    symbol: str,
    h_electron_affinity: float = DEFAULT_H_ELECTRON_AFFINITY,
) -> AtomParameters:
    """Electronegativity X and hardness J for an element, per EQeq's formulas."""
    if symbol == "H":
        # X = 0.5 * (hI1 + hI0); J = hI1 - hI0
        x = 0.5 * (_HI1 + h_electron_affinity)
        j = _HI1 - h_electron_affinity
        return AtomParameters(electronegativity=x, hardness=j)

    if symbol not in _ELEMENT_INDEX:
        raise ValueError(f"element {symbol!r} is outside EQeq's parameter table")
    z = _ELEMENT_INDEX[symbol]
    ip = _ionization_table()[z]
    cc = _charge_centers().get(symbol, 0)
    if cc + 1 >= len(ip):
        raise ValueError(f"missing ionization data for {symbol} at charge center {cc}")

    # X = 0.5 * (IP[cc+1] + IP[cc]); J = IP[cc+1] - IP[cc]; X -= cc * J
    x = 0.5 * (ip[cc + 1] + ip[cc])
    j = ip[cc + 1] - ip[cc]
    x -= cc * j
    return AtomParameters(electronegativity=x, hardness=j)
