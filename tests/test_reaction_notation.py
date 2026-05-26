import pytest

from reaction_kinematics import Reaction


def test_notation_equivalent_to_four_args():
    """3H(p,n)3He notation gives the same reaction as four explicit masses."""
    rxn_notation = Reaction("3H(p,n)3He")
    rxn_explicit = Reaction("p", "3H", "n", "3He")
    assert rxn_notation.q_value == pytest.approx(rxn_explicit.q_value)


def test_notation_with_spaces():
    """Spaces around the comma and parentheses are accepted."""
    rxn = Reaction("3H( p , n )3He")
    assert rxn.q_value == pytest.approx(Reaction("p", "3H", "n", "3He").q_value)


def test_notation_elastic_scattering():
    """Elastic scattering notation a(12C,a)12C works and has Q=0."""
    rxn_notation = Reaction("12C(a,a)12C")
    rxn_explicit = Reaction("a", "12C", "a", "12C")
    assert rxn_notation.q_value == pytest.approx(rxn_explicit.q_value)
    assert rxn_notation.q_value == pytest.approx(0.0, abs=1e-6)


def test_notation_kinematics_match():
    """Kinematic results from notation and explicit construction are identical."""
    rxn_notation = Reaction("3H(p,n)3He")
    rxn_explicit = Reaction("p", "3H", "n", "3He")
    table_n = rxn_notation.kinematics_table_at_beam_energy(1.2)
    table_e = rxn_explicit.kinematics_table_at_beam_energy(1.2)
    import numpy as np
    for key in table_n:
        assert np.allclose(table_n[key], table_e[key])


def test_notation_invalid_raises():
    """Malformed notation string (has parenthesis but wrong format) raises ValueError."""
    with pytest.raises(ValueError, match="Invalid reaction notation"):
        Reaction("3H(p)3He")


def test_notation_mixed_args_raises():
    """Mixing notation string with extra mass args raises ValueError."""
    with pytest.raises(ValueError, match="Cannot mix"):
        Reaction("3H(p,n)3He", "extra")
