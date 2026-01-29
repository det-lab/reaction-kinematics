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


rxn = TwoBody("alpha", "12C", "alpha", '12C' ,4 , mass_unit="MeV")
      
kinematics_arrays = rxn.compute_arrays()


# exact

# interpolated single value
print(rxn.at_value("theta3", 25 * math.pi/180, y_names="e3"))

# interpolated multiple outputs
print(rxn.at_value(
    "theta3",
    25 * math.pi/180,
    y_names=["e3", "v3", "p3"]
))

# interpolated full state
print(rxn.at_value("theta_cm", 0.8))


import matplotlib.pyplot as plt

# assuming rxn is already created

theta3 = kinematics_arrays["theta3"]   # radians
e3 = kinematics_arrays["e3"]           # energy (same units as input)

plt.figure()
plt.plot(e3, theta3)
plt.xlabel(r"$\theta_3$ (rad)")
plt.ylabel(r"$E_3$")
plt.title("Ejectile Energy vs Angle")
plt.grid(True)
plt.show()
