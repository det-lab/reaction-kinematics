#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 26 13:02:45 2025

@author: joey
"""

import sys
from pathlib import Path

# Absolute path to THIS file
HERE = Path(__file__).resolve()

# Go up to Nuclear/
NUCLEAR = HERE.parents[1]

# Path to your package's src directory
SRC = NUCLEAR / "reaction-kinematics" / "src"

sys.path.insert(0, str(SRC))

print("Using SRC:", SRC)
print("sys.path[0]:", sys.path[0])

from reaction_kinematics.inputs import MassInput, EnergyValue
from reaction_kinematics.units import EnergyUnit
from reaction_kinematics.reaction_kinematics import TwoBody


import reaction_kinematics
print(reaction_kinematics.__file__)

# ----------------------------
# Simple physical sanity test
# ----------------------------

# Elastic scattering: p + 12C → p + 12C
p = MassInput("p")
C12 = MassInput("12C")

ek = EnergyValue(5.0, EnergyUnit.MeV)

rxn = TwoBody(
    p.mass,
    C12.mass,
    p.mass,
    C12.mass,
    ek.value
)

print("\n--- Reaction Test: p + 12C ---")
print("s =", rxn.s)
print("Emax ejectile (MeV) =", rxn.emax3)
print("Emin ejectile (MeV) =", rxn.emin3)

if rxn.theta3max is not None:
    print("Theta_max (deg) =", rxn.theta3max * 180 / math.pi)
else:
    print("Theta_max: none (full angular range)")
