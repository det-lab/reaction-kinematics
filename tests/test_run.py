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
    data = rxn.kinematics_table_at_beam_energy(1.2)

    # Basic interpolation at a fixed variable value
    result = rxn.kinematics_at_beam_energy_and_angle(1.2, "theta_cm", 46)
    assert isinstance(result, dict)
    assert "theta_cm" in result  # theta_cm key name is unchanged
    assert all(isinstance(v, float) for values in result.values() for v in values)

    # Interpolation near the edge of the theta4 range
    theta4_arr = data["theta4_lab"]
    edge_angle = max(theta4_arr) - 1e-12
    for y in ("energy3_lab", "velocity3_lab", "momentum3_lab"):
        vals = rxn.kinematics_at_beam_energy_and_angle(1.2, "theta4_lab", edge_angle)[y]
        assert all(isinstance(v, float) for v in vals)

    # Plot energy vs angle without errors
    theta4 = data["theta4_lab"]
    e3 = data["energy3_lab"]
    plt.figure()
    plt.plot(theta4, e3)
    plt.xlabel(r"$\theta_4$ (deg)")
    plt.ylabel(r"$E_3$")
    plt.title("Ejectile Energy vs Angle")
    plt.grid(True)
    plt.close()
