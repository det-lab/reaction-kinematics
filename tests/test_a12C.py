#!/usr/bin/env python3
"""
Created on Fri Feb 27 21:08:26 2026

@author: joey
"""

from pathlib import Path

import numpy as np
import pandas as pd

from reaction_kinematics import Reaction


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

    rxn = Reaction("a", "12C", "a", "12C")

    calc = []

    for _, row in df.iterrows():
        theta3 = float(np.deg2rad(row["theta3 (degrees)"]))
        r = rxn.at_value("theta3", theta3, ek=4.0)

        # at_value can return multiple theta4 branches, so choose the one
        # that matches the reference table value for this row
        theta4_vals_deg = np.rad2deg(np.array(r["theta4"]))
        theta4_best = theta4_vals_deg[np.argmin(np.abs(theta4_vals_deg - row["theta4 (degrees)"]))]

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
    theta4_mask = ~np.isclose(df["theta3 (degrees)"], 0.0, atol=1e-12)

    assert np.allclose(calc["E3"], df["E3 (MeV)"], rtol=1e-3)
    assert np.allclose(calc["E4"], df["E4 (MeV)"], rtol=1e-3)
    assert np.allclose(calc["theta3cm"], df["theta3cm (degrees)"], rtol=1e-3)
    assert np.allclose(
        calc.loc[theta4_mask, "theta4"], df.loc[theta4_mask, "theta4 (degrees)"], rtol=1e-3
    )
    assert np.allclose(calc["v3"], df["v3 (fraction of c)"], rtol=1e-3)
    assert np.allclose(calc["v4"], df["v4 (fraction of c)"], rtol=1e-3)
