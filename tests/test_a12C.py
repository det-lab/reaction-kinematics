#!/usr/bin/env python3
"""
Created on Fri Feb 27 21:08:26 2026

@author: joey
"""

from pathlib import Path

import numpy as np
import pandas as pd

from reaction_kinematics import TwoBody


def test_alpha12C_reference_table():
    """
    Validate TwoBody results against published alpha + 12C elastic scattering data.
    Source:
    https://skisickness.com/math-physics/kinematics/
    Retrieved March 5 2026
    """

    # Load reference data
    data_file = Path(__file__).parent / "data" / "alpha_12C_alpha_12C_reference_data.csv"

    df = pd.read_csv(data_file, skiprows=6)

    rxn = TwoBody("a", "12C", "a", "12C", 4.0, mass_unit="MeV")

    calc = []

    for _, row in df.iterrows():
        theta3 = np.deg2rad(row["theta3"])
        r = rxn.at_value("theta3", theta3)

        theta4_vals_deg = np.rad2deg(np.array(r["theta4"]))
        theta4_best = theta4_vals_deg[np.argmin(np.abs(theta4_vals_deg - row["theta4"]))]

        calc.append(
            {
                "E3": r["e3"][0],
                "E4": r["e4"][0],
                "theta3cm": np.rad2deg(r["theta_cm"][0]),
                "theta4": theta4_best,
                "v3": r["v3"][0],
                "v4": r["v4"][0],
            }
        )

    calc = pd.DataFrame(calc)

    # vectorized comparisons
    theta4_mask = ~np.isclose(df["theta3"], 0.0, atol=1e-12)

    assert np.allclose(calc["E3"], df["E3"], rtol=1e-3)
    assert np.allclose(calc["E4"], df["E4"], rtol=1e-3)
    assert np.allclose(calc["theta3cm"], df["theta3cm"], rtol=1e-3)
    assert np.allclose(calc.loc[theta4_mask, "theta4"], df.loc[theta4_mask, "theta4"], rtol=1e-3)
    assert np.allclose(calc["v3"], df["v3"], rtol=1e-3)
    assert np.allclose(calc["v4"], df["v4"], rtol=1e-3)
