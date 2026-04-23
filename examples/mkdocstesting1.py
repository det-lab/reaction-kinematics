#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 10:20:04 2026

@author: joey
"""
from reaction_kinematics import TwoBody
import matplotlib.pyplot as plt

rxn = TwoBody("p", "3H", "n", "3He", 1.2)
#Proton + Tritium Reaction

data = rxn.compute_arrays()

plt.plot(data["theta4"], data["e3"])
plt.xlabel("Recoil Angle θ₄ (rad)")
plt.ylabel("Ejectile Energy E₃ (MeV)")
plt.title("E₃ vs θ₄")
plt.grid(True)
plt.show()