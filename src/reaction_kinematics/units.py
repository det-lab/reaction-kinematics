"""
Unit definitions and conversions
"""

import math
from enum import Enum


class EnergyUnit(Enum):
    keV = 1e-3
    MeV = 1.0
    GeV = 1e3
    TeV = 1e6

    @classmethod
    def from_any(cls, unit):
        if isinstance(unit, cls):
            return unit
        if isinstance(unit, str):
            try:
                return cls[unit]
            except KeyError:
                raise ValueError(f"Unknown energy unit '{unit}'")
        raise TypeError(f"Invalid unit type: {type(unit)}")


class AngleUnit(Enum):
    rad = 1.0
    deg = math.pi / 180.0
    mrad = 1e-3
