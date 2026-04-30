"""
User input containers
"""

from .constants import AMU
from .mass import get_mass
from .units import EnergyUnit


class MassInput:
    def __init__(self, value, unit="ael"):
        if unit == "ael":
            self.mass = get_mass(value)
        elif unit == "MeV":
            self.mass = float(value)
        elif unit == "amu":
            self.mass = float(value) * AMU
        else:
            raise ValueError(f"Invalid mass unit: {unit}")


class EnergyValue:
    def __init__(self, value, unit=EnergyUnit.MeV):
        unit = EnergyUnit.from_any(unit)
        self.value = float(value) * unit.value
