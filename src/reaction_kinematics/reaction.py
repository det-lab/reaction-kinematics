"""
Relativistic two-body reaction kinematics
"""

import math
from collections.abc import Iterable

import numpy as np
import numpy.typing as npt

from reaction_kinematics.inputs import MassInput
from reaction_kinematics.units import AngleUnit, EnergyUnit

MassArg = str | int | float | MassInput


def _parse_mass(m: MassArg, unit: str | EnergyUnit | None = None) -> float:
    if isinstance(m, MassInput):
        return m.mass
    if isinstance(m, str):
        return MassInput(m).mass
    if isinstance(m, (int, float)):
        if unit is None:
            raise ValueError("Numeric mass provided without mass_unit")
        if isinstance(unit, str):
            unit = EnergyUnit[unit]
        return m * unit.value
    raise TypeError(f"Unsupported mass input type: {type(m)}")


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
        Number of CM angle grid points (default 500). Changing this invalidates
        the cache.

    Examples
    --------
    >>> rxn = Reaction("p", "3H", "n", "3He")
    >>> rxn.q_value
    -0.763...
    >>> data = rxn.compute_arrays(ek=1.2)
    >>> result = rxn.at_value("theta3", 0.5, ek=1.2)
    >>> branches = rxn.kinematic_curve(np.deg2rad(30), np.linspace(1.0, 5.0, 200))
    """

    def __init__(
        self,
        m1: MassArg,
        m2: MassArg,
        m3: MassArg,
        m4: MassArg,
        *,
        mass_unit: str | EnergyUnit | None = None,
    ) -> None:
        self._m1 = _parse_mass(m1, mass_unit)
        self._m2 = _parse_mass(m2, mass_unit)
        self._m3 = _parse_mass(m3, mass_unit)
        self._m4 = _parse_mass(m4, mass_unit)
        self._ncoscm: int = 500
        self._cached_ek: float | None = None
        self._table: dict[str, list[float]] | None = None
        # per-energy kinematic state, populated by _bind
        self._nogo: bool = False
        self._pcmp: float | None = None
        self._thesinh: float | None = None
        self._thecosh: float | None = None
        self._e03: float | None = None
        self._e04: float | None = None
        # lab-frame energy extrema
        self._emax3: float | None = None
        self._emin3: float | None = None
        self._emax4: float | None = None
        self._emin4: float | None = None
        # max lab angles and associated quantities (None = no forward maximum)
        self._theta3max: float | None = None
        self._theta4max: float | None = None
        self._e3atmaxang: float | None = None
        self._e4atmaxang: float | None = None
        self._cmcos3max: float | None = None
        self._cmcos4max: float | None = None

    @property
    def _ncoscm(self) -> int:
        return self._ncoscm

    @_ncoscm.setter
    def _ncoscm(self, val: int) -> None:
        self._ncoscm = val
        self._cached_ek = None
        self._table = None

    @property
    def q_value(self) -> float:
        """Q-value of the reaction in MeV: Q = m1 + m2 - m3 - m4."""
        return self._m1 + self._m2 - self._m3 - self._m4

    def _bind(self, ek_mev: float) -> None:
        """Compute and cache kinematic state for the given beam energy."""
        if ek_mev == self._cached_ek:
            return
        self._cached_ek = ek_mev
        self._table = None
        self._compute(ek_mev)

    def _compute(self, ek_mev: float) -> None:
        self._nogo = False
        self._pcmp = self._thesinh = self._thecosh = self._e03 = self._e04 = None

        # Mandelstam s
        s = (self._m1 + self._m2) ** 2 + 2.0 * self._m2 * ek_mev
        if s <= 0.0:
            self._nogo = True
            return

        # initial CM momentum
        pcm2 = (s - self._m1**2 - self._m2**2) ** 2 - 4.0 * self._m1**2 * self._m2**2
        if pcm2 < 0:
            self._nogo = True
            return
        pcm = math.sqrt(pcm2 / (4.0 * s))

        # CM rapidity → boost parameters
        acmratio = (math.sqrt(self._m2**2 + pcm**2) + pcm) / self._m2
        cmrap = math.log(acmratio)
        self._thesinh = math.sinh(cmrap)
        self._thecosh = math.cosh(cmrap)

        # final-state CM momentum
        pcmp2 = (s - self._m3**2 - self._m4**2) ** 2 - 4.0 * self._m3**2 * self._m4**2
        if pcmp2 < 0:
            self._nogo = True
            return
        self._pcmp = math.sqrt(pcmp2 / (4.0 * s))

        # CM total energies of outgoing particles
        self._e03 = math.sqrt(self._pcmp**2 + self._m3**2)
        self._e04 = math.sqrt(self._pcmp**2 + self._m4**2)

        # lab-frame energy extrema
        self._emax3 = self._e03 * self._thecosh + self._pcmp * self._thesinh - self._m3
        self._emin3 = self._e03 * self._thecosh - self._pcmp * self._thesinh - self._m3
        self._emax4 = self._e04 * self._thecosh + self._pcmp * self._thesinh - self._m4
        self._emin4 = self._e04 * self._thecosh - self._pcmp * self._thesinh - self._m4

        # reset max-angle quantities to sentinel values
        self._theta3max = None
        self._theta4max = None
        self._e3atmaxang = None
        self._e4atmaxang = None
        self._cmcos3max = None
        self._cmcos4max = None

        # max ejectile lab angle (only exists when pcmp < m3 * sinh(rapidity))
        thetatest3: float | None = None
        if self._m3 > 0.0:
            thetatest3 = self._pcmp / (self._m3 * self._thesinh)
            if thetatest3 < 1.0:
                self._theta3max = math.asin(thetatest3)
                patmax = (self._e03 * math.cos(self._theta3max) * self._thesinh) / (
                    1.0 + thetatest3**2 * self._thesinh**2
                )
                eatmax = math.sqrt(patmax**2 + self._m3**2)
                self._e3atmaxang = eatmax - self._m3
                self._cmcos3max = (eatmax - self._e03 * self._thecosh) / (
                    self._pcmp * self._thesinh
                )

        # elastic case: forward-angle symmetry forces theta3max = 90°
        if (self._m1 + self._m2) == (self._m3 + self._m4) and thetatest3 is not None:
            if abs(thetatest3 - 1.0) < 1e-3:
                self._theta3max = math.pi / 2.0
                self._cmcos3max = -1.0
                self._e3atmaxang = (
                    self._e03 * self._thecosh
                    + self._cmcos3max * self._pcmp * self._thesinh
                    - self._m3
                )

        # max recoil lab angle
        thetatest4: float | None = None
        if self._m4 > 0.0:
            thetatest4 = self._pcmp / (self._m4 * self._thesinh)
            if thetatest4 < 1.0:
                self._theta4max = math.asin(thetatest4)
                patmax = (self._e04 * math.cos(self._theta4max) * self._thesinh) / (
                    1.0 + thetatest4**2 * self._thesinh**2
                )
                eatmax = math.sqrt(patmax**2 + self._m4**2)
                self._e4atmaxang = eatmax - self._m4
                self._cmcos4max = (eatmax - self._e04 * self._thecosh) / (
                    self._pcmp * self._thesinh
                )

        # elastic case: theta4max = 90°
        if (self._m1 + self._m2) == (self._m3 + self._m4) and thetatest4 is not None:
            if abs(thetatest4 - 1.0) < 1e-3:
                self._theta4max = math.pi / 2.0
                self._cmcos4max = 1.0
                self._e4atmaxang = (
                    self._e04 * self._thecosh
                    - self._cmcos4max * self._pcmp * self._thesinh
                    - self._m4
                )

    def _kinematics_at_coscm(self, coscm: float) -> dict[str, float]:
        if (
            self._pcmp is None
            or self._thecosh is None
            or self._e03 is None
            or self._thesinh is None
            or self._e04 is None
        ):
            raise ValueError("Kinematic state not computed — call _bind first")

        sincm = math.sqrt(max(0.0, 1.0 - coscm**2))

        ppar3 = self._pcmp * self._thecosh * coscm + self._e03 * self._thesinh
        pperp3 = self._pcmp * sincm
        ptot3 = math.hypot(ppar3, pperp3)

        ppar4 = -self._pcmp * self._thecosh * coscm + self._e04 * self._thesinh
        pperp4 = self._pcmp * sincm
        ptot4 = math.hypot(ppar4, pperp4)

        e3 = self._e03 * self._thecosh + coscm * self._pcmp * self._thesinh - self._m3
        e4 = self._e04 * self._thecosh - coscm * self._pcmp * self._thesinh - self._m4

        return {
            "coscm": coscm,
            "theta_cm": math.acos(coscm),
            "theta3": math.acos(ppar3 / ptot3) if ptot3 > 0 else 0.0,
            "theta4": math.acos(ppar4 / ptot4) if ptot4 > 0 else 0.0,
            "e3": e3,
            "e4": e4,
            "v3": ptot3 / (e3 + self._m3),
            "v4": ptot4 / (e4 + self._m4),
            "p3": ptot3,
            "p4": ptot4,
        }

    def _build_table(self) -> None:
        """Build the interpolation table over the full CM angle grid."""
        keys = ["coscm", "theta_cm", "theta3", "theta4", "e3", "e4", "v3", "v4", "p3", "p4"]
        table: dict[str, list[float]] = {k: [] for k in keys}
        for i in range(-self._ncoscm, self._ncoscm + 1):
            row = self._kinematics_at_coscm(i / self._ncoscm)
            for k in keys:
                table[k].append(row[k])
        self._table = table

    def compute_arrays(
        self, ek: float, *, energy_unit: EnergyUnit = EnergyUnit.MeV
    ) -> dict[str, npt.NDArray[np.float64]]:
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
        dict[str, np.ndarray]
            Keys: ``"coscm"``, ``"theta_cm"``, ``"theta3"``, ``"theta4"``,
            ``"e3"``, ``"e4"``, ``"v3"``, ``"v4"``.

        Raises
        ------
        ValueError
            If the reaction is kinematically forbidden at this energy.
        """
        ek_mev = _parse_energy(ek, energy_unit)
        self._bind(ek_mev)
        if self._nogo:
            raise ValueError(f"Reaction kinematically forbidden at ek={ek_mev} MeV")
        keys = ["coscm", "theta_cm", "theta3", "theta4", "e3", "e4", "v3", "v4"]
        rows = [
            self._kinematics_at_coscm(i / self._ncoscm)
            for i in range(-self._ncoscm, self._ncoscm + 1)
        ]
        return {k: np.array([row[k] for row in rows]) for k in keys}

    def at_value(
        self,
        x_name: str,
        x: float,
        *,
        ek: float,
        energy_unit: EnergyUnit = EnergyUnit.MeV,
        y_names: str | list[str] | None = None,
        duplicate_tol: float = 1e-6,
    ) -> dict[str, list[float]]:
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
        self._bind(ek_mev)

        if self._table is None:
            self._build_table()
        assert self._table is not None

        xs = self._table[x_name]

        if isinstance(y_names, str):
            y_names = [y_names]
        if y_names is None:
            y_names = list(self._table.keys())

        solutions = []

        exact_idx = np.where(np.isclose(xs, x, atol=1e-12))[0]
        if len(exact_idx) > 0:
            for i in exact_idx:
                solutions.append({k: self._table[k][i] for k in y_names})
        else:
            found = False
            for i in range(len(xs) - 1):
                x0, x1 = xs[i], xs[i + 1]
                if (x0 - x) * (x1 - x) <= 0 and x0 != x1:
                    found = True
                    t = (x - x0) / (x1 - x0)
                    solutions.append(
                        {
                            k: self._table[k][i] + t * (self._table[k][i + 1] - self._table[k][i])
                            for k in y_names
                        }
                    )
            if not found:
                raise ValueError(f"{x_name}={x} outside physical range")

        e_key = "e3" if "e3" in y_names else y_names[0]
        unique: list = []
        for sol in solutions:
            if not any(abs(sol[e_key] - u[e_key]) < duplicate_tol for u in unique):
                unique.append(sol)
        unique.sort(key=lambda s: s[e_key], reverse=True)

        return {k: [s[k] for s in unique] for k in y_names}

    def kinematic_curve(
        self,
        theta: float,
        ek_array: Iterable[float],
        *,
        angle_unit: AngleUnit = AngleUnit.rad,
        energy_unit: EnergyUnit = EnergyUnit.MeV,
    ) -> list[dict[str, npt.NDArray[np.float64]]]:
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
            try:
                row = self.at_value("theta3", theta_rad, ek=ek_mev, y_names=keys)
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
