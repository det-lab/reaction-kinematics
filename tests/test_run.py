#!/usr/bin/env python3
"""
Created on Fri Dec 26 13:02:45 2025

@author: joey
"""

# Smoke-test example and plotting in a pytest function
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from reaction_kinematics import TwoBody


def test_run_and_plot_example() -> None:
    """
    Smoke-test the TwoBody kinematics example and ensure plotting works.
    """
    # Initialize reaction and compute arrays
    rxn = TwoBody("a", "12C", "a", "12C", 1.2, mass_unit="MeV")
    data = rxn.compute_arrays()

    # Basic at_value interpolation
    result = rxn.at_value("theta_cm", 0.8)
    assert isinstance(result, dict)
    assert "theta_cm" in result
    assert all(isinstance(v, float) for values in result.values() for v in values)

    # Ensure theta4max is defined and interpolation at its edge does not error
    assert rxn.theta4max is not None
    edge_angle = rxn.theta4max - 9.76e-8
    for y in ("e3", "v3", "p3"):
        vals = rxn.at_value("theta4", edge_angle, y_names=y)[y]
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
