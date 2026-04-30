#!/usr/bin/env python3
"""
Created on Fri Dec 26 13:02:45 2025

@author: joey
"""

# Smoke-test example and plotting in a pytest function
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from reaction_kinematics import Reaction


def test_run_and_plot_example() -> None:
    """
    Smoke-test the Reaction kinematics example and ensure plotting works.
    """
    rxn = Reaction("a", "12C", "a", "12C")
    data = rxn.compute_arrays(1.2)

    # Basic at_value interpolation
    result = rxn.at_value("theta_cm", 0.8, ek=1.2)
    assert isinstance(result, dict)
    assert "theta_cm" in result
    assert all(isinstance(v, float) for values in result.values() for v in values)

    # Interpolation near the edge of the theta4 range
    theta4_arr = data["theta4"]
    edge_angle = max(theta4_arr) - 1e-12
    for y in ("e3", "v3", "p3"):
        vals = rxn.at_value("theta4", edge_angle, ek=1.2, y_names=y)[y]
        assert all(isinstance(v, float) for v in vals)

    # Plot energy vs angle without errors
    theta4 = data["theta4"]
    e3 = data["e3"]
    plt.figure()
    plt.plot(theta4, e3)
    plt.xlabel(r"$\theta_4$ (rad)")
    plt.ylabel(r"$E_3$")
    plt.title("Ejectile Energy vs Angle")
    plt.grid(True)
    plt.close()
