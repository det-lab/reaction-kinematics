"""
Relativistic two-body reaction kinematics
"""

import numpy as np

from reaction_kinematics.two_body import TwoBody as _TwoBody
from reaction_kinematics.two_body import _parse_mass
from reaction_kinematics.units import AngleUnit, EnergyUnit


def _parse_energy(ek: float, energy_unit: EnergyUnit) -> float:
    if isinstance(energy_unit, str):
        energy_unit = EnergyUnit[energy_unit]
    return float(ek) * energy_unit.value


class Reaction:
    """
    Defines a two-body nuclear reaction: projectile + target → ejectile + recoil.

    All energy-dependent methods accept an ``ek`` parameter (beam kinetic energy,
    MeV by default). Internal computation is cached per energy so repeated calls
    at the same energy are efficient.

    Parameters
    ----------
    m1, m2, m3, m4 : str, MassInput, or float
        Masses of projectile, target, ejectile, and recoil.
        Strings like ``"p"``, ``"12C"``, ``"alpha"`` are looked up in the mass table.
        Floats require ``mass_unit``.
    mass_unit : str or EnergyUnit, optional
        Unit for numeric masses (e.g. ``"MeV"``, ``"keV"``).

    Attributes
    ----------
    ncoscm : int
        Number of CM angle grid points used by :meth:`compute_arrays` and
        :meth:`at_value` (default 500). Changing this invalidates the cache.

    Examples
    --------
    >>> rxn = Reaction("p", "3H", "n", "3He")
    >>> rxn.q_value
    -0.763...
    >>> data = rxn.compute_arrays(ek=1.2)
    >>> result = rxn.at_value("theta3", 0.5, ek=1.2)
    >>> branches = rxn.kinematic_curve(np.deg2rad(30), np.linspace(1.0, 5.0, 200))
    """

    def __init__(self, m1, m2, m3, m4, *, mass_unit=None):
        self.m1 = _parse_mass(m1, mass_unit)
        self.m2 = _parse_mass(m2, mass_unit)
        self.m3 = _parse_mass(m3, mass_unit)
        self.m4 = _parse_mass(m4, mass_unit)
        self._ncoscm: int = 500
        self._cached_ek: float | None = None
        self._bound: _TwoBody | None = None

    @property
    def ncoscm(self) -> int:
        return self._ncoscm

    @ncoscm.setter
    def ncoscm(self, val: int) -> None:
        self._ncoscm = val
        self._cached_ek = None

    @property
    def q_value(self) -> float:
        """Q-value of the reaction in MeV: Q = m1 + m2 - m3 - m4."""
        return self.m1 + self.m2 - self.m3 - self.m4

    def _at_energy(self, ek_mev: float) -> _TwoBody:
        if ek_mev != self._cached_ek:
            self._bound = _TwoBody(self.m1, self.m2, self.m3, self.m4, ek_mev, mass_unit="MeV")
            self._bound.ncoscm = self._ncoscm
            self._cached_ek = ek_mev
        assert self._bound is not None
        return self._bound

    def compute_arrays(self, ek: float, *, energy_unit: EnergyUnit = EnergyUnit.MeV) -> dict:
        """
        Compute full kinematics over a CM angle grid.

        Parameters
        ----------
        ek : float
            Beam kinetic energy.
        energy_unit : EnergyUnit, optional
            Unit of ``ek`` (default MeV).

        Returns
        -------
        dict of str -> list
            Keys: ``"coscm"``, ``"theta_cm"``, ``"theta3"``, ``"theta4"``,
            ``"e3"``, ``"e4"``, ``"v3"``, ``"v4"``.

        Raises
        ------
        ValueError
            If the reaction is kinematically forbidden at this energy.
        """
        ek_mev = _parse_energy(ek, energy_unit)
        bound = self._at_energy(ek_mev)
        if bound.nogo:
            raise ValueError(f"Reaction kinematically forbidden at ek={ek_mev} MeV")
        return bound.compute_arrays()

    def at_value(
        self,
        x_name: str,
        x: float,
        *,
        ek: float,
        energy_unit: EnergyUnit = EnergyUnit.MeV,
        y_names=None,
        duplicate_tol: float = 1e-6,
    ) -> dict:
        """
        Interpolate kinematic quantities at a fixed independent variable.

        Always returns lists to handle multi-valued cases (e.g. two ejectile
        energies at the same lab angle).

        Parameters
        ----------
        x_name : str
            Independent variable name, e.g. ``"theta3"``, ``"theta_cm"``,
            ``"coscm"``, ``"theta4"``.
        x : float
            Value to evaluate at (radians for angles).
        ek : float
            Beam kinetic energy.
        energy_unit : EnergyUnit, optional
            Unit of ``ek`` (default MeV).
        y_names : str or list of str, optional
            Dependent variables to return. ``None`` returns all.
        duplicate_tol : float, optional
            Tolerance for merging near-duplicate solutions (default 1e-6).

        Returns
        -------
        dict of str -> list
            Each value is a list of solutions, sorted descending by ``e3``.

        Raises
        ------
        ValueError
            If ``x`` is outside the physical range.

        Examples
        --------
        >>> rxn = Reaction("p", "3H", "n", "3He")
        >>> rxn.at_value("theta3", 0.5, ek=1.2, y_names=["e3", "v3"])
        {'e3': [...], 'v3': [...]}
        """
        ek_mev = _parse_energy(ek, energy_unit)
        return self._at_energy(ek_mev).at_value(
            x_name, x, y_names=y_names, duplicate_tol=duplicate_tol
        )

    def kinematic_curve(
        self,
        theta: float,
        ek_array,
        *,
        angle_unit: AngleUnit = AngleUnit.rad,
        energy_unit: EnergyUnit = EnergyUnit.MeV,
    ) -> list[dict]:
        """
        Compute ejectile kinematics at a fixed lab angle over a range of beam energies.

        Returns two branches (high- and low-energy) as a list of two dicts. Each
        dict contains arrays indexed by beam energy, with ``NaN`` where that branch
        does not exist. Branch 0 is always the higher-energy solution.

        Parameters
        ----------
        theta : float
            Fixed lab angle of the ejectile.
        ek_array : array-like
            Beam energies to sweep.
        angle_unit : AngleUnit, optional
            Unit of ``theta`` (default radians).
        energy_unit : EnergyUnit, optional
            Unit of ``ek_array`` values (default MeV).

        Returns
        -------
        list of two dicts, each with keys:
            ``"ek"``, ``"e3"``, ``"e4"``, ``"theta4"``, ``"v3"``, ``"v4"``.

        Examples
        --------
        >>> rxn = Reaction("p", "3H", "n", "3He")
        >>> branches = rxn.kinematic_curve(np.deg2rad(30), np.linspace(1.0, 5.0, 200))
        >>> for b in branches:
        ...     plt.plot(b["ek"], b["e3"])
        """
        if isinstance(angle_unit, str):
            angle_unit = AngleUnit[angle_unit]
        theta_rad = theta * angle_unit.value

        keys = ["e3", "e4", "theta4", "v3", "v4"]
        branches = [
            {"ek": [], **{k: [] for k in keys}},
            {"ek": [], **{k: [] for k in keys}},
        ]

        for ek in ek_array:
            ek_mev = _parse_energy(ek, energy_unit)
            bound = self._at_energy(ek_mev)

            try:
                row = bound.at_value("theta3", theta_rad, y_names=keys)
            except ValueError:
                solutions = []
            else:
                n = len(row["e3"])
                solutions = [{k: row[k][i] for k in keys} for i in range(n)]

            for i, branch in enumerate(branches):
                branch["ek"].append(ek_mev)
                sol = solutions[i] if i < len(solutions) else None
                for k in keys:
                    branch[k].append(sol[k] if sol is not None else float("nan"))

        return [{k: np.array(v) for k, v in branch.items()} for branch in branches]
