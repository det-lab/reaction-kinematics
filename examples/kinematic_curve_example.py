#!/usr/bin/env python3
"""
Kinematic curve for 3H(p,n)3He at a fixed lab angle of 30 degrees.

Plots neutron energy vs proton beam energy. Near threshold the reaction
is double-valued; both branches are plotted.
"""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from reaction_kinematics import Reaction

ek_array = np.linspace(1.0, 5.0, 500)
branches = Reaction("p", "3H", "n", "3He").kinematic_curve(np.deg2rad(30), ek_array)

fig, ax = plt.subplots(figsize=(7, 5))

for branch in branches:
    ax.plot(branch["ek"], branch["e3"])

ax.set_xlabel("Proton beam energy $E_p$ (MeV)")
ax.set_ylabel("Neutron energy $E_n$ (MeV)")
ax.set_title(r"$^3$H(p,n)$^3$He at 30°")
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig("kinematic_curve_plot.png", dpi=150)
print("Saved kinematic_curve_plot.png")
