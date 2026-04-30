#!/usr/bin/env python3
"""
Created on Thu Mar 12 10:51:29 2026

@author: joey
"""

from pathlib import Path

import numpy as np
import pandas as pd

from reaction_kinematics import Reaction


def test_q_value_p3H_n3He():
    """
    Validate Q-value for 3H(p,n)3He against NNDC Q-value calculator.
    Reference: https://www.nndc.bnl.gov/qcalc/
    Q = -763.75498 keV = -0.76375498 MeV
    """
    expected_mev = -763.75498e-3  # keV → MeV

    rxn = Reaction("p", "3H", "n", "3He")
    assert np.isclose(rxn.q_value, expected_mev, rtol=1e-5)


def test_p3H_n3He_reference_table():
    """
    Validate TwoBody results against published 3H(p,n)3He reference data.

    The reference file may contain rows that correspond to multiple valid
    kinematic branches, so for multi-valued outputs we select the branch
    closest to the reference value in that row.
    """

    data_file = Path(__file__).parent / "data" / "p_3H_n_3He_reference_data.csv"
    df = pd.read_csv(data_file, skiprows=6)

    rxn = Reaction("p", "3H", "n", "3He")

    calc = []

    for _, row in df.iterrows():
        r = rxn.at_value("coscm", float(row["costheta3cm"]), ek=1.2)  # pyright: ignore[reportArgumentType]

        e3_vals = np.array(r["e3"])
        e4_vals = np.array(r["e4"])
        theta_cm_vals_deg = np.rad2deg(np.array(r["theta_cm"]))
        theta4_vals_deg = np.rad2deg(np.array(r["theta4"]))
        v3_vals = np.array(r["v3"])
        v4_vals = np.array(r["v4"])

        # Some inputs can return multiple valid branches.
        # Pick the solution closest to the published reference row.
        e3_best = e3_vals[np.argmin(np.abs(e3_vals - row["E3 (MeV)"]))]
        e4_best = e4_vals[np.argmin(np.abs(e4_vals - row["E4 (MeV)"]))]
        theta_cm_best = theta_cm_vals_deg[
            np.argmin(np.abs(theta_cm_vals_deg - row["theta3cm (degrees)"]))
        ]
        theta4_best = theta4_vals_deg[np.argmin(np.abs(theta4_vals_deg - row["theta4 (degrees)"]))]
        v3_best = v3_vals[np.argmin(np.abs(v3_vals - row["v3 (fraction of c)"]))]
        v4_best = v4_vals[np.argmin(np.abs(v4_vals - row["v4 (fraction of c)"]))]

        calc.append(
            {
                "E3 (MeV)": e3_best,
                "E4 (MeV)": e4_best,
                "theta3cm (degrees)": theta_cm_best,
                "theta4 (degrees)": theta4_best,
                "v3 (fraction of c)": v3_best,
                "v4 (fraction of c)": v4_best,
            }
        )

    calc = pd.DataFrame(calc)

    assert np.allclose(calc["E3 (MeV)"], df["E3 (MeV)"], rtol=1e-3)
    assert np.allclose(calc["E4 (MeV)"], df["E4 (MeV)"], rtol=1e-3)
    assert np.allclose(calc["theta3cm (degrees)"], df["theta3cm (degrees)"], rtol=1e-3)
    assert np.allclose(calc["theta4 (degrees)"], df["theta4 (degrees)"], rtol=1e-3)
    assert np.allclose(calc["v3 (fraction of c)"], df["v3 (fraction of c)"], rtol=1e-3)
    assert np.allclose(calc["v4 (fraction of c)"], df["v4 (fraction of c)"], rtol=1e-3)
