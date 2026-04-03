#!/usr/bin/env python3
"""
Created on Thu Mar 19 09:47:14 2026

@author: joey
"""

import pytest

from reaction_kinematics import TwoBody


def test_p3H_n3He_below_threshold_raises():
    """
    3H(p,n)3He should not produce valid kinematics at 0.5 MeV.
    """
    with pytest.raises(ValueError):
        TwoBody("p", "3H", "n", "3He", 0.5, mass_unit="MeV").compute_arrays()
