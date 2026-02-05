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


rxn = TwoBody("p", "3H", "n", '3He' ,1.2, mass_unit="MeV")
      
kinematics_arrays = rxn.compute_arrays()
# interpolated full state
print(rxn.at_value("theta_cm", 0.8))
eps = 9.76e-8   # small number

angle = rxn.theta4max - eps

# exact

# interpolated single value
print(rxn.at_value("theta4", angle, y_names="e3"))
print('----------')

# interpolated multiple outputs
print(rxn.at_value(
    "theta4",
    angle,
    y_names=["e3", "v3", "p3"], duplicate_tol = 1e-3
))
# interpolated full state



import matplotlib.pyplot as plt

# assuming rxn is already created

theta4 = kinematics_arrays["theta4"]   # radians
e3 = kinematics_arrays["e3"]           # energy (same units as input)

plt.figure()
plt.plot(theta4, e3 )
plt.ylabel(r"$\theta_4$ (rad)")
plt.xlabel(r"$E_3$")
plt.title("Ejectile Energy vs Angle")
plt.grid(True)
plt.show()
