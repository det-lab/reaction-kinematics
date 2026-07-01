"""
Tests for kinematics_curve_at_angle against digitized data from Krane "Introductory Nuclear Physics".
Data files encode the reaction and fixed lab angle in the filename.
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from reaction_kinematics import Reaction

DATA_DIR = Path(__file__).parent / "data"

KRANE_FILES = sorted(DATA_DIR.glob("3H(p,n)3He * (Krane) *.csv"))

# Maximum acceptable deviation from digitized Krane plot values (MeV).
# Reading off a printed plot introduces ~60 keV of digitization uncertainty.
KRANE_ATOL_MEV = 0.06


def parse_angle(path: Path) -> float:
    """Extract lab angle in degrees from filename."""
    match = re.search(r"(\d+) degrees", path.name)
    if not match:
        raise ValueError(f"Could not parse angle from filename: {path.name}")
    return float(match.group(1))


def test_kinematic_curve_return_structure():
    """kinematics_curve_at_angle returns a list of two dicts each with the expected keys."""
    beam_energy_array = np.linspace(1.5, 5.0, 20)
    branches = Reaction("p", "3H", "n", "3He").kinematics_curve_at_angle(beam_energy_array, 30)

    assert len(branches) == 2
    expected_keys = {
        "beam_energy_lab",
        "energy3_lab",
        "energy4_lab",
        "theta4_lab",
        "velocity3_lab",
        "velocity4_lab",
        "momentum3_lab",
        "momentum4_lab",
    }
    for b in branches:
        assert set(b.keys()) == expected_keys
        for v in b.values():
            assert len(v) == len(beam_energy_array)


def test_kinematic_curve_single_valued():
    """3H(p,n)3He at 30° is single-valued over this energy range: branch 0 fully
    populated, branch 1 all NaN."""
    beam_energy_array = np.linspace(1.5, 5.0, 50)
    branches = Reaction("p", "3H", "n", "3He").kinematics_curve_at_angle(beam_energy_array, 30)

    assert not np.any(np.isnan(branches[0]["energy3_lab"]))
    assert np.all(np.isnan(branches[1]["energy3_lab"]))


def test_kinematic_curve_two_valued():
    """12C(p,p)12C in inverse kinematics at 3° has a two-valued region: both
    branches should be populated for most of the energy range."""
    # Two-valued regime for 12C(p,p)12C at 3° starts well below 50 MeV 12C beam energy.
    beam_energy_array = np.linspace(50.0, 200.0, 100)
    branches = Reaction("12C", "p", "12C", "p").kinematics_curve_at_angle(beam_energy_array, 3)

    # Branch 0 should be fully populated
    assert not np.any(np.isnan(branches[0]["energy3_lab"]))
    # Branch 1 should have values for more than 90% of the range
    assert np.sum(~np.isnan(branches[1]["energy3_lab"])) > 0.9 * len(beam_energy_array)
    # Branch 0 always higher energy than branch 1 where both exist
    both = ~np.isnan(branches[1]["energy3_lab"])
    assert np.all(branches[0]["energy3_lab"][both] > branches[1]["energy3_lab"][both])


@pytest.mark.parametrize("data_file", KRANE_FILES, ids=lambda p: p.name)
def test_kinematic_curve_vs_krane(data_file):
    """
    Validate kinematics_curve_at_angle E3 output against digitized Krane data.

    Near threshold some angles are double-valued, so each reference point is
    matched to the nearest computed E3 from either branch. Absolute tolerance
    of 60 keV accounts for digitization error from reading values off a printed
    plot; rtol is avoided since E3 values near threshold can be very small.
    """
    theta3_deg = parse_angle(data_file)

    df = pd.read_csv(data_file, header=None, names=["beam_energy_lab", "e3_ref"])

    branches = Reaction("p", "3H", "n", "3He").kinematics_curve_at_angle(
        df["beam_energy_lab"].to_numpy(), theta3_deg
    )

    e3_b0 = branches[0]["energy3_lab"]
    e3_b1 = branches[1]["energy3_lab"]

    for i, e3_ref in enumerate(df["e3_ref"]):
        # Skip points where both branches are NaN — either below threshold
        # or at the resolution limit near theta3max.
        if np.isnan(e3_b0[i]) and np.isnan(e3_b1[i]):
            continue

        candidates = [v for v in [e3_b0[i], e3_b1[i]] if not np.isnan(v)]
        best = min(candidates, key=lambda v: abs(v - e3_ref))
        assert np.isclose(best, e3_ref, atol=KRANE_ATOL_MEV, rtol=0.0), (
            f"ek={df['beam_energy_lab'][i]:.4f}: computed {best:.4f}, reference {e3_ref:.4f} "
            f"(angle={theta3_deg} deg, file={data_file.name})"
        )
