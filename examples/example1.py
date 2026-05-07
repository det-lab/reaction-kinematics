#!/usr/bin/env python3
"""
Smoke-test the Reaction API and basic plotting.
"""

# Smoke-test example and plotting in a pytest function
import matplotlib
print(matplotlib.get_backend())


import matplotlib.pyplot as plt

from reaction_kinematics import Reaction


"""
Smoke-test the Reaction kinematics example and ensure plotting works.
"""
# Initialize reaction and compute arrays
rxn = Reaction("alpha", "12C", "alpha", "12C")
data = rxn.kinematics_table_at_beam_energy(4.0)

# Basic interpolation at a fixed variable value
result = rxn.kinematics_at_beam_energy_and_variable(4.0, "theta_cm", 0.8)
assert isinstance(result, dict)
assert "theta_cm" in result
assert all(isinstance(v, float) for values in result.values() for v in values)

# Ensure theta4max is defined and interpolation at its edge does not error
# assert rxn.theta4max is not None
# edge_angle = rxn.theta4max - 9.76e-8
# for y in ("e3", "v3", "p3"):
   # vals = rxn.kinematics_at_beam_energy_and_variable(1.2, "theta4", edge_angle, return_variables=y)[y]
   # assert all(isinstance(v, float) for v in vals)

# Plot energy vs angle without errors
theta4 = data["theta4"]
e3 = data["e3"]
plt.figure()
plt.plot(theta4, e3)
plt.xlabel(r"$\theta_4$ (rad)")
plt.ylabel(r"$E_3$")
plt.title("Ejectile Energy vs Angle")
plt.grid(True)
plt.show()