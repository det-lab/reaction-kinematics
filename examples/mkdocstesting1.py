#!/usr/bin/env python3
"""
Plotting example: Ejectile Energy vs Recoil Angle for p + 3H → n + 3He.
"""
import matplotlib.pyplot as plt

from reaction_kinematics import Reaction

rxn = Reaction("p", "3H", "n", "3He")
data = rxn.kinematics_table_at_beam_energy(1.2)

plt.plot(data["theta4"], data["e3"])
plt.xlabel("Recoil Angle θ₄ (rad)")
plt.ylabel("Ejectile Energy E₃ (MeV)")
plt.title("E₃ vs θ₄")
plt.grid(True)
plt.show()
