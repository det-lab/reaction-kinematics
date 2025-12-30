#!/usr/bin/env python3
"""
Created on Fri Dec 26 13:02:45 2025

@author: joey
"""

import math

import reaction_kinematics
from reaction_kinematics.inputs import EnergyValue, MassInput
from reaction_kinematics.reaction_kinematics import TwoBody
from reaction_kinematics.units import EnergyUnit

print(reaction_kinematics.__file__)

# ----------------------------
# Simple physical sanity test
# ----------------------------

# Elastic scattering: p + 12C → p + 12C
p = MassInput("p")
C12 = MassInput("12C")

ek = EnergyValue(5.0, EnergyUnit.MeV)

rxn = TwoBody(p.mass, C12.mass, p.mass, C12.mass, ek.value)

print("\n--- Reaction Test: p + 12C ---")
print("s =", rxn.s)
print("Emax ejectile (MeV) =", rxn.emax3)
print("Emin ejectile (MeV) =", rxn.emin3)

if rxn.theta3max is not None:
    print("Theta_max (deg) =", rxn.theta3max * 180 / math.pi)
else:
    print("Theta_max: none (full angular range)")
