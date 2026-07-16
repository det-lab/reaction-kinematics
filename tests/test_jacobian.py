"""
Tests for the analytic dOmega_lab/dOmega_cm Jacobian (jacobian3_lab, jacobian4_lab)
against independent numerical differentiation of the lab/cm angle relationship.

Numerical differentiation is only accurate with a much finer cm-angle grid than
the library's default (1001 points): np.gradient's finite-difference error is
otherwise as large as a few percent in steep regions of theta_lab(theta_cm).
"""

import numpy as np

from reaction_kinematics import Reaction
from reaction_kinematics.units import AngleUnit


def _numeric_jacobian(cos_cm: np.ndarray, theta_lab: np.ndarray) -> np.ndarray:
    return np.gradient(np.cos(theta_lab), cos_cm)


def test_jacobian3_matches_numeric_derivative():
    rxn = Reaction("p", "11B", "a", "8Be")
    rxn.n_cm_grid_points = 20001
    table = rxn.kinematics_table_at_beam_energy(5.4, angle_unit=AngleUnit.rad)

    numeric = _numeric_jacobian(table["cos_theta_cm"], table["theta3_lab"])
    analytic = table["jacobian3_lab"]

    # Exclude points near the grid edges, where np.gradient's one-sided
    # finite difference is a poor approximation of the true derivative.
    interior = slice(400, -400)
    assert np.allclose(analytic[interior], numeric[interior], atol=1e-3)


def test_jacobian4_matches_numeric_derivative():
    rxn = Reaction("p", "11B", "a", "8Be")
    rxn.n_cm_grid_points = 20001
    table = rxn.kinematics_table_at_beam_energy(5.4, angle_unit=AngleUnit.rad)

    numeric = _numeric_jacobian(table["cos_theta_cm"], table["theta4_lab"])
    analytic = table["jacobian4_lab"]

    interior = slice(400, -400)
    assert np.allclose(analytic[interior], numeric[interior], atol=1e-3)


def test_jacobians_match_numeric_derivative_elastic_case():
    """3H(p,n)3He at a beam energy well above threshold, single-valued kinematics."""
    rxn = Reaction("p", "3H", "n", "3He")
    rxn.n_cm_grid_points = 20001
    table = rxn.kinematics_table_at_beam_energy(5.0, angle_unit=AngleUnit.rad)

    interior = slice(400, -400)
    numeric3 = _numeric_jacobian(table["cos_theta_cm"], table["theta3_lab"])
    numeric4 = _numeric_jacobian(table["cos_theta_cm"], table["theta4_lab"])

    assert np.allclose(table["jacobian3_lab"][interior], numeric3[interior], atol=1e-3)
    assert np.allclose(table["jacobian4_lab"][interior], numeric4[interior], atol=1e-3)
