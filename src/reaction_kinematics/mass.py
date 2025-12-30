"""
Mass lookup and isotope parsing
"""

import re

from .constants import AMU, EMASS
from .masstable import MTAB


def parse_isotope(isotope: str):
    """
    Parse strings like '12C', '4He', 'p', 'd', 't'
    Returns (A, element)
    """
    isotope = isotope.strip().lower()

    aliases = {
        "p": (1, "h"),
        "d": (2, "h"),
        "t": (3, "h"),
        "alpha": (4, "he"),
    }

    if isotope in aliases:
        return aliases[isotope]

    match = re.match(r"(\d+)([a-z]+)", isotope)
    if not match:
        raise ValueError(f"Invalid isotope string: {isotope}")

    return int(match.group(1)), match.group(2)


def get_mass(isotope: str) -> float:
    """
    Return nuclear mass in MeV
    """
    A, el = parse_isotope(isotope)

    try:
        entry = MTAB[(A, el)]
    except KeyError as err:
        raise KeyError(f"Isotope not found: {A}{el}") from err

    mass_excess = entry["mexcess"] * 1e-3  # keV → MeV
    return A * AMU + mass_excess - entry["Z"] * EMASS
