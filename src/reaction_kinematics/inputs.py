"""
User input containers
"""

from .mass import get_mass
from .units import AngleUnit, EnergyUnit


class MassInput:
    def __init__(self, value, unit="ael"):
        if unit == "ael":
            self.mass = get_mass(value)
        elif unit == "MeV":
            self.mass = float(value)
        elif unit == "amu":
            self.mass = float(value) * 931.494104
        else:
            raise ValueError(f"Invalid mass unit: {unit}")


class EnergyValue:
    def __init__(self, value, unit: EnergyUnit):
        self.value = float(value) * unit.value


class ReactionInput:
    def __init__(
        self,
        projectile,
        target,
        ejectile,
        recoil,
        kinetic_energy: EnergyValue,
        ex_ejectile: EnergyValue | None = None,
        ex_recoil: EnergyValue | None = None,
        angle_unit: AngleUnit = AngleUnit.deg,
        energy_unit: EnergyUnit = EnergyUnit.MeV,
    ):
        if ex_ejectile is None:
            ex_ejectile = EnergyValue(0.0, EnergyUnit.MeV)
        if ex_recoil is None:
            ex_recoil = EnergyValue(0.0, EnergyUnit.MeV)

        self.m1 = projectile.mass
        self.m2 = target.mass
        self.m3 = ejectile.mass
        self.m4 = recoil.mass

        self.ek = kinetic_energy.value
        self.ex3 = ex_ejectile.value if ex_ejectile else 0.0
        self.ex4 = ex_recoil.value if ex_recoil else 0.0

        self.angle_unit = angle_unit
        self.energy_unit = energy_unit
