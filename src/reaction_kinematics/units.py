"""
Unit definitions and conversions
"""

import math
from enum import Enum
from typing import TypeVar

T = TypeVar("T", bound="_UnitEnum")


class _UnitEnum(Enum):
    """Enum mixin adding a friendly ``from_any(enum_member_or_name)`` constructor."""

    @classmethod
    def from_any(cls: type[T], unit: "T | str") -> T:
        if isinstance(unit, cls):
            return unit
        if isinstance(unit, str):
            try:
                return cls[unit]
            except KeyError as err:
                valid = ", ".join(m.name for m in cls)
                raise ValueError(
                    f"Unknown {cls.__name__} '{unit}' (expected one of: {valid})"
                ) from err
        raise TypeError(f"Invalid unit type: {type(unit)}")


class EnergyUnit(_UnitEnum):
    keV = 1e-3
    MeV = 1.0
    GeV = 1e3
    TeV = 1e6


class AngleUnit(_UnitEnum):
    rad = 1.0
    deg = math.pi / 180.0
    mrad = 1e-3
