"""
Tests that energy_unit governs every energy- and momentum-valued output
(not just the beam_energy input), mirroring how angle_unit already governs
every angle-valued output. Also covers output_units().
"""

import numpy as np
import pytest

from reaction_kinematics import Reaction
from reaction_kinematics.inputs import EnergyValue, MassInput
from reaction_kinematics.units import AngleUnit, EnergyUnit, ureg

ENERGY_KEYS = ["energy3_lab", "energy4_lab", "momentum3_lab", "momentum4_lab"]
DIMENSIONLESS_KEYS = [
    "cos_theta_cm",
    "velocity3_lab",
    "velocity4_lab",
    "jacobian3_lab",
    "jacobian4_lab",
]


def test_kinematics_table_energy_unit_scales_outputs():
    rxn = Reaction("p", "11B", "a", "8Be")

    # Same physical beam energy expressed in each unit.
    t_mev = rxn.kinematics_table_at_beam_energy(5.4, energy_unit=EnergyUnit.MeV)
    t_kev = rxn.kinematics_table_at_beam_energy(5400.0, energy_unit=EnergyUnit.keV)

    for key in ENERGY_KEYS:
        assert np.allclose(t_kev[key], 1000.0 * t_mev[key])

    for key in DIMENSIONLESS_KEYS + ["theta3_lab", "theta4_lab", "theta_cm"]:
        assert np.allclose(t_kev[key], t_mev[key])


def test_kinematics_at_beam_energy_and_angle_energy_unit_scales_outputs():
    rxn = Reaction("p", "11B", "a", "8Be")

    r_mev = rxn.kinematics_at_beam_energy_and_angle(
        5.4, "theta3_lab", 16.07, energy_unit=EnergyUnit.MeV
    )
    r_kev = rxn.kinematics_at_beam_energy_and_angle(
        5400.0, "theta3_lab", 16.07, energy_unit=EnergyUnit.keV
    )

    for key in ENERGY_KEYS:
        assert np.allclose(r_kev[key], [1000.0 * v for v in r_mev[key]])


def test_kinematics_curve_at_angle_energy_unit_scales_outputs():
    rxn = Reaction("p", "11B", "a", "8Be")

    branches_mev = rxn.kinematics_curve_at_angle(
        np.linspace(4.0, 6.0, 5), 16.07, energy_unit=EnergyUnit.MeV
    )
    branches_kev = rxn.kinematics_curve_at_angle(
        np.linspace(4000.0, 6000.0, 5), 16.07, energy_unit=EnergyUnit.keV
    )

    assert np.allclose(
        branches_kev[0]["beam_energy_lab"], 1000.0 * branches_mev[0]["beam_energy_lab"]
    )
    for key in ENERGY_KEYS:
        assert np.allclose(branches_kev[0][key], 1000.0 * branches_mev[0][key])


def test_duplicate_tol_is_interpreted_in_energy_unit():
    """12C(p,p)12C at 3 degrees, inverse kinematics: two branches ~22 MeV apart."""
    rxn = Reaction("12C", "p", "12C", "p")

    # A 25 MeV tolerance merges the two branches into one solution.
    merged = rxn.kinematics_at_beam_energy_and_angle(100.0, "theta3_lab", 3.0, duplicate_tol=25)
    assert len(merged["energy3_lab"]) == 1

    # The same tolerance value, but as 25 keV (=0.025 MeV), should not merge
    # the (still ~22 MeV apart) branches.
    unmerged = rxn.kinematics_at_beam_energy_and_angle(
        100_000.0, "theta3_lab", 3.0, duplicate_tol=25, energy_unit=EnergyUnit.keV
    )
    assert len(unmerged["energy3_lab"]) == 2


def test_output_units_reflects_requested_units():
    default_units = Reaction.output_units()
    assert default_units["theta3_lab"] == ureg.Unit("deg")
    assert default_units["energy3_lab"] == ureg.Unit("MeV")
    assert default_units["momentum3_lab"] == ureg.Unit("MeV/c")
    assert default_units["cos_theta_cm"] == ureg.Unit("dimensionless")
    assert default_units["velocity3_lab"] == ureg.Unit("dimensionless")

    custom_units = Reaction.output_units(angle_unit="rad", energy_unit="keV")
    assert custom_units["theta3_lab"] == ureg.Unit("rad")
    assert custom_units["energy3_lab"] == ureg.Unit("keV")
    assert custom_units["momentum3_lab"] == ureg.Unit("keV/c")


def test_kinematics_result_units_match_output_units_and_never_drift():
    """The whole point of bundling .units with the data: they're built in the
    same call from the same angle_unit/energy_unit, so they can't disagree."""
    rxn = Reaction("p", "11B", "a", "8Be")
    expected = Reaction.output_units(angle_unit=AngleUnit.rad, energy_unit=EnergyUnit.keV)

    table = rxn.kinematics_table_at_beam_energy(
        5400.0, angle_unit=AngleUnit.rad, energy_unit=EnergyUnit.keV
    )
    for key in table:
        assert table.units[key] == expected[key]

    result = rxn.kinematics_at_beam_energy_and_angle(
        5400.0, "theta3_lab", 0.28, angle_unit=AngleUnit.rad, energy_unit=EnergyUnit.keV
    )
    for key in result:
        assert result.units[key] == expected[key]

    branches = rxn.kinematics_curve_at_angle(
        np.linspace(4000.0, 6000.0, 5),
        0.28,
        angle_unit=AngleUnit.rad,
        energy_unit=EnergyUnit.keV,
    )
    for key in branches[0]:
        assert branches[0].units[key] == expected[key]


@pytest.mark.parametrize(
    "method_keys",
    [
        {
            "cos_theta_cm",
            "theta_cm",
            "theta3_lab",
            "theta4_lab",
            "energy3_lab",
            "energy4_lab",
            "velocity3_lab",
            "velocity4_lab",
            "momentum3_lab",
            "momentum4_lab",
            "jacobian3_lab",
            "jacobian4_lab",
        },
        {
            "beam_energy_lab",
            "energy3_lab",
            "energy4_lab",
            "theta4_lab",
            "velocity3_lab",
            "velocity4_lab",
            "momentum3_lab",
            "momentum4_lab",
            "jacobian3_lab",
            "jacobian4_lab",
        },
    ],
)
def test_output_units_covers_every_result_key(method_keys):
    units = Reaction.output_units()
    assert method_keys <= units.keys()


def test_outside_physical_range_error_shows_original_value_and_unit():
    """12C(p,p)12C at 3 degrees has a ~4.8 degree max lab angle; 10 degrees
    is unreachable. The error must echo what was typed (10 deg), not the
    internally-converted radians value."""
    rxn = Reaction("12C", "p", "12C", "p")
    with pytest.raises(ValueError, match=r"theta3_lab=10 deg outside physical range"):
        rxn.kinematics_at_beam_energy_and_angle(100.0, "theta3_lab", 10, angle_unit=AngleUnit.deg)


def test_outside_physical_range_error_respects_angle_unit():
    rxn = Reaction("12C", "p", "12C", "p")
    with pytest.raises(ValueError, match=r"theta3_lab=1\.5 rad outside physical range"):
        rxn.kinematics_at_beam_energy_and_angle(100.0, "theta3_lab", 1.5, angle_unit=AngleUnit.rad)


def test_cos_theta_cm_outside_range_error_has_no_spurious_unit():
    rxn = Reaction("p", "11B", "a", "8Be")
    with pytest.raises(ValueError, match=r"^cos_theta_cm=5\.0 outside physical range$"):
        rxn.kinematics_at_beam_energy_and_angle(5.4, "cos_theta_cm", 5.0)


@pytest.mark.parametrize(
    "call,expected_substr",
    [
        (lambda: Reaction(float("nan"), "11B", "a", "8Be", mass_unit="keV"), "nan keV"),
        (lambda: MassInput(float("nan"), unit="amu"), "nan amu"),
        (lambda: EnergyValue(float("nan"), unit="GeV"), "nan GeV"),
    ],
)
def test_non_finite_input_errors_show_the_unit_it_was_given_in(call, expected_substr):
    with pytest.raises(ValueError, match=r"is not a finite number") as exc_info:
        call()
    assert expected_substr in str(exc_info.value)


def test_non_finite_beam_energy_error_shows_the_unit_it_was_given_in():
    rxn = Reaction("p", "11B", "a", "8Be")
    with pytest.raises(ValueError, match=r"nan keV is not a finite number"):
        rxn.kinematics_at_beam_energy_and_angle(
            float("nan"), "theta3_lab", 16.07, energy_unit=EnergyUnit.keV
        )


def test_non_finite_angle_value_error_shows_the_unit_it_was_given_in():
    rxn = Reaction("p", "11B", "a", "8Be")
    with pytest.raises(ValueError, match=r"nan rad is not a finite number"):
        rxn.kinematics_at_beam_energy_and_angle(
            5.4, "theta3_lab", float("nan"), angle_unit=AngleUnit.rad
        )
