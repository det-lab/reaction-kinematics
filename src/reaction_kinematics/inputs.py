"""
User input containers
"""

from .constants import AMU
from .mass import get_mass
from .units import EnergyUnit


class MassInput:
    mass: float

    def __init__(self, value: str | float, unit: str = "ael") -> None:
        if unit == "ael":
            self.mass = get_mass(value)  # type: ignore[arg-type]
        elif unit == "MeV":
            self.mass = float(value)
        elif unit == "amu":
            self.mass = float(value) * AMU
        else:
            raise ValueError(f"Invalid mass unit: {unit}")


class EnergyValue:
    value: float

    def __init__(self, value: float, unit: EnergyUnit | str = EnergyUnit.MeV) -> None:
        unit = EnergyUnit.from_any(unit)
        self.value = float(value) * unit.value
