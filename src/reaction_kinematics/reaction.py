"""
Relativistic two-body reaction kinematics
"""

import math
import re
from collections.abc import Iterable, Iterator, Mapping
from typing import TypeVar

import numpy as np
import numpy.typing as npt
import pint

from reaction_kinematics.inputs import MassInput
from reaction_kinematics.units import AngleUnit, EnergyUnit, ureg

_V = TypeVar("_V")

MassArg = str | int | float | MassInput

# Keys in a kinematics result dict that are angle-valued (scaled by angle_unit)
# or energy-valued (scaled by energy_unit) on output. Everything else
# (cos_theta_cm, velocityN_lab, jacobianN_lab) is dimensionless.
_ANGLE_KEYS = frozenset({"theta_cm", "theta3_lab", "theta4_lab"})
_MOMENTUM_KEYS = frozenset({"momentum3_lab", "momentum4_lab"})
_ENERGY_KEYS = frozenset({"energy3_lab", "energy4_lab", "beam_energy_lab"}) | _MOMENTUM_KEYS

# Every key that can appear in a dict returned by kinematics_table_at_beam_energy,
# kinematics_at_beam_energy_and_angle, or kinematics_curve_at_angle.
_ALL_RESULT_KEYS = (
    _ANGLE_KEYS
    | _ENERGY_KEYS
    | frozenset(
        {"cos_theta_cm", "velocity3_lab", "velocity4_lab", "jacobian3_lab", "jacobian4_lab"}
    )
)


def _result_unit(key: str, angle_unit: AngleUnit, energy_unit: EnergyUnit) -> pint.Unit:
    """The pint Unit for a kinematics result key, e.g. Unit('MeV/c') for momentum."""
    if key in _ANGLE_KEYS:
        return ureg.Unit(angle_unit.name)
    if key in _MOMENTUM_KEYS:
        return ureg.Unit(f"{energy_unit.name}/c")
    if key in _ENERGY_KEYS:
        return ureg.Unit(energy_unit.name)
    return ureg.Unit("dimensionless")


def _result_units(
    keys: Iterable[str], angle_unit: AngleUnit, energy_unit: EnergyUnit
) -> dict[str, pint.Unit]:
    return {k: _result_unit(k, angle_unit, energy_unit) for k in keys}


class KinematicsResult(Mapping[str, _V]):
    """
    A read-only, dict-like kinematics result: ``result["theta3_lab"]`` works
    exactly like a plain dict of arrays/lists, plus a ``.units`` dict (str ->
    ``pint.Unit``) giving the unit of each key.

    Values are plain numbers/arrays, not pint.Quantity — they work unmodified
    with numpy, matplotlib, pandas, etc. Check ``.units[key]`` when you need to
    know or convert the unit; wrap a value yourself (``ureg.Quantity(result[key],
    result.units[key])``) if you want pint-native arithmetic for a specific step.
    """

    def __init__(self, data: dict[str, _V], units: dict[str, pint.Unit]) -> None:
        self._data = data
        self.units = units

    def __getitem__(self, key: str) -> _V:
        return self._data[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._data)

    def __len__(self) -> int:
        return len(self._data)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self._data!r}, units={self.units!r})"


def _parse_mass(m: MassArg, unit: str | EnergyUnit | None = None) -> float:
    if isinstance(m, MassInput):
        return m.mass
    if isinstance(m, str):
        return MassInput(m).mass
    if isinstance(m, (int, float)):
        if unit is None:
            raise ValueError("Numeric mass provided without mass_unit")
        unit = EnergyUnit.from_any(unit)
        mass = m * unit.value
        if not math.isfinite(mass):
            raise ValueError(f"mass={m} {unit.name} is not a finite number")
        return mass
    raise TypeError(f"Unsupported mass input type: {type(m)}")


def _parse_reaction_notation(notation: str) -> tuple[str, str, str, str]:
    """Parse 'target(beam,ejectile)recoil' → (beam, target, ejectile, recoil)."""
    match = re.fullmatch(r"\s*(\S+)\s*\(\s*(\S+)\s*,\s*(\S+)\s*\)\s*(\S+)\s*", notation)
    if not match:
        raise ValueError(
            f"Invalid reaction notation {notation!r}. "
            "Expected 'target(beam,ejectile)recoil', e.g. '3H(p,n)3He'."
        )
    target, beam, ejectile, recoil = match.groups()
    return beam, target, ejectile, recoil


def _parse_energy(ek: float, energy_unit: EnergyUnit) -> float:
    energy_unit = EnergyUnit.from_any(energy_unit)
    ek_mev = float(ek) * energy_unit.value
    if not math.isfinite(ek_mev):
        raise ValueError(f"beam_energy={ek} {energy_unit.name} is not a finite number")
    return ek_mev


class Reaction:
    """
    Defines a two-body nuclear reaction: projectile + target → ejectile + recoil.

    All energy-dependent methods accept a ``beam_energy`` parameter (beam kinetic
    energy, MeV by default) and an ``energy_unit`` keyword that governs both that
    input and every energy- and momentum-valued output; likewise ``angle_value``/
    ``angle_unit`` govern every angle-valued input and output.

    The three ``kinematics_*`` methods below return a ``KinematicsResult``: a
    read-only, dict-like object (``result["theta3_lab"]`` works like a plain
    dict of arrays/lists) with a ``.units`` dict (str -> ``pint.Unit``) giving
    the unit of each key. ``Reaction.output_units(angle_unit=..., energy_unit=...)``
    gives that same ``{key: pint.Unit}`` mapping without needing to run a
    computation first. See the ``Returns`` section of each method below for
    per-key detail. Internal computation is cached per energy so repeated
    calls at the same energy are efficient.

    Parameters
    ----------
    mass1 : str, MassInput, or float
        Projectile mass, or a full reaction string in ``"target(beam,ejectile)recoil"``
        notation (e.g. ``"3H(p,n)3He"``). If notation is given, ``mass2``–``mass4``
        must be omitted.
    mass2, mass3, mass4 : str, MassInput, or float, optional
        Target, ejectile, and recoil masses. Required when ``mass1`` is not a
        reaction notation string. Strings like ``"p"``, ``"12C"``, ``"alpha"``
        are looked up in the mass table (no ``mass_unit`` needed — the table is
        already in the right units). Floats require ``mass_unit``.
    mass_unit : str or EnergyUnit, optional
        Unit for numeric masses (e.g. ``"MeV"``, ``"keV"``). Only used for
        ``float``/``int`` mass arguments; ignored for isotope-string arguments.

    Attributes
    ----------
    n_cm_grid_points : int
        Total number of CM angle grid points (default 1001). Changing this
        invalidates the cache.

    Examples
    --------
    >>> rxn = Reaction("p", "3H", "n", "3He")
    >>> rxn = Reaction("3H(p,n)3He")  # equivalent
    >>> rxn.q_value
    -0.763...
    >>> data = rxn.kinematics_table_at_beam_energy(1.2)
    >>> result = rxn.kinematics_at_beam_energy_and_angle(1.2, "theta3_lab", 30)
    >>> branches = rxn.kinematics_curve_at_angle(np.linspace(1.0, 5.0, 200), 30)
    """

    def __init__(
        self,
        mass1: MassArg,
        mass2: MassArg | None = None,
        mass3: MassArg | None = None,
        mass4: MassArg | None = None,
        *,
        mass_unit: str | EnergyUnit | None = None,
    ) -> None:
        if isinstance(mass1, str) and "(" in mass1:
            if any(m is not None for m in (mass2, mass3, mass4)):
                raise ValueError(
                    "Cannot mix reaction notation string with separate mass arguments."
                )
            _m1, _m2, _m3, _m4 = _parse_reaction_notation(mass1)
        else:
            if mass2 is None or mass3 is None or mass4 is None:
                raise ValueError("Provide either a reaction notation string or all four masses.")
            _m1, _m2, _m3, _m4 = mass1, mass2, mass3, mass4
        self._m1 = _parse_mass(_m1, mass_unit)
        self._m2 = _parse_mass(_m2, mass_unit)
        self._m3 = _parse_mass(_m3, mass_unit)
        self._m4 = _parse_mass(_m4, mass_unit)
        self.__n_cm_grid_points: int = 1001
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
    def n_cm_grid_points(self) -> int:
        return self.__n_cm_grid_points

    @n_cm_grid_points.setter
    def n_cm_grid_points(self, val: int) -> None:
        if not isinstance(val, int):
            raise TypeError(f"n_cm_grid_points must be an int, got {type(val).__name__}")
        if val < 2:
            raise ValueError(f"n_cm_grid_points={val} must be at least 2")
        self.__n_cm_grid_points = val
        self._cached_ek = None
        self._table = None

    @property
    def q_value(self) -> float:
        """Q-value of the reaction in MeV: Q = (beam + target) - (ejectile + recoil)."""
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

        # dOmega_lab/dOmega_cm for each ejectile, i.e. the factor that converts a
        # lab-frame differential cross section to the cm-frame one:
        # dsigma/dOmega_cm = dsigma/dOmega_lab * jacobianN_lab
        #
        # Both are d(cos_theta_lab)/d(coscm), coscm being cos_theta_cm of particle 3.
        # ppar4's sign convention (see above) flips the sign of the thecosh term
        # relative to particle 3's formula.
        if ptot3 > 0:
            jacobian3_lab = (self._pcmp**2 / ptot3**2) * (
                coscm * ppar3 / ptot3 + sincm * pperp3 * self._thecosh / ptot3
            )
        else:
            jacobian3_lab = math.nan
        if ptot4 > 0:
            jacobian4_lab = (self._pcmp**2 / ptot4**2) * (
                coscm * ppar4 / ptot4 - sincm * pperp4 * self._thecosh / ptot4
            )
        else:
            jacobian4_lab = math.nan

        return {
            "cos_theta_cm": coscm,
            "theta_cm": math.acos(coscm),
            "theta3_lab": math.acos(ppar3 / ptot3) if ptot3 > 0 else 0.0,
            "theta4_lab": math.acos(ppar4 / ptot4) if ptot4 > 0 else 0.0,
            "energy3_lab": e3,
            "energy4_lab": e4,
            "velocity3_lab": ptot3 / (e3 + self._m3),
            "velocity4_lab": ptot4 / (e4 + self._m4),
            "momentum3_lab": ptot3,
            "momentum4_lab": ptot4,
            "jacobian3_lab": jacobian3_lab,
            "jacobian4_lab": jacobian4_lab,
        }

    def _build_table(self) -> None:
        """Build the interpolation table over the full CM angle grid."""
        keys = [
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
        ]
        table: dict[str, list[float]] = {k: [] for k in keys}
        for coscm in np.linspace(-1.0, 1.0, self.n_cm_grid_points):
            row = self._kinematics_at_coscm(coscm)
            for k in keys:
                table[k].append(row[k])
        self._table = table

    def kinematics_table_at_beam_energy(
        self,
        beam_energy: float,
        *,
        angle_unit: AngleUnit = AngleUnit.deg,
        energy_unit: EnergyUnit = EnergyUnit.MeV,
    ) -> KinematicsResult[npt.NDArray[np.float64]]:
        """
        Compute full kinematics over a CM angle grid.

        Values are plain arrays (work unmodified with numpy/matplotlib/etc.).
        The returned ``KinematicsResult`` also carries a ``.units`` dict (str ->
        ``pint.Unit``) giving the unit of each key.

        Parameters
        ----------
        beam_energy : float
            Beam kinetic energy, in ``energy_unit``.
        angle_unit : AngleUnit, optional
            Unit for angle outputs ``theta_cm``, ``theta3_lab``, ``theta4_lab``
            (default degrees). Does not affect ``jacobianN_lab``, which is a
            dimensionless ratio of solid angles regardless of angle unit.
        energy_unit : EnergyUnit, optional
            Unit of ``beam_energy`` (default MeV). Also governs the unit of every
            energy- and momentum-valued output (``energy3_lab``, ``energy4_lab``,
            ``momentum3_lab``, ``momentum4_lab``), exactly as ``angle_unit`` governs
            every angle-valued output.

        Returns
        -------
        KinematicsResult[np.ndarray]
            ``"cos_theta_cm"`` : dimensionless, in [-1, 1].
            ``"theta_cm"``, ``"theta3_lab"``, ``"theta4_lab"`` : ``angle_unit``
            (default degrees).
            ``"energy3_lab"``, ``"energy4_lab"`` : ``energy_unit`` (default MeV).
            ``"velocity3_lab"``, ``"velocity4_lab"`` : dimensionless, as a fraction
            of the speed of light ``c``.
            ``"momentum3_lab"``, ``"momentum4_lab"`` : ``energy_unit``/``c``
            (natural units, default MeV/``c``).
            ``"jacobian3_lab"``, ``"jacobian4_lab"`` : dimensionless.
            ``jacobianN_lab`` is ``dOmegaN_lab/dOmega_cm``, i.e. the factor that
            converts a lab-frame differential cross section to the cm-frame one:
            ``dsigma/dOmega_cm = dsigma/dOmegaN_lab * jacobianN_lab``.

            ``result.units`` carries this same mapping as ``{key: pint.Unit}``.
            ``Reaction.output_units(angle_unit=..., energy_unit=...)`` gives the
            same thing without needing to run a computation first.

        Raises
        ------
        ValueError
            If the reaction is kinematically forbidden at this energy.
        """
        angle_unit = AngleUnit.from_any(angle_unit)
        energy_unit = EnergyUnit.from_any(energy_unit)
        ek_mev = _parse_energy(beam_energy, energy_unit)
        self._bind(ek_mev)
        if self._nogo:
            raise ValueError(f"Reaction kinematically forbidden at beam_energy={ek_mev} MeV")
        keys = [
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
        ]
        rows = [
            self._kinematics_at_coscm(coscm)
            for coscm in np.linspace(-1.0, 1.0, self.n_cm_grid_points)
        ]
        result = {k: np.array([row[k] for row in rows]) for k in keys}
        for k in keys:
            if k in _ANGLE_KEYS:
                result[k] = result[k] / angle_unit.value
            elif k in _ENERGY_KEYS:
                result[k] = result[k] / energy_unit.value
        return KinematicsResult(result, _result_units(keys, angle_unit, energy_unit))

    def kinematics_at_beam_energy_and_angle(
        self,
        beam_energy: float,
        angle_name: str,
        angle_value: float,
        *,
        angle_unit: AngleUnit = AngleUnit.deg,
        energy_unit: EnergyUnit = EnergyUnit.MeV,
        duplicate_tol: float = 1e-6,
    ) -> KinematicsResult[list[float]]:
        """
        Interpolate kinematic quantities at a fixed beam energy and kinematic variable value.

        Always returns lists to handle multi-valued cases (e.g. two ejectile
        energies at the same lab angle). Values are plain lists (work unmodified
        with numpy/pandas/etc.). The returned ``KinematicsResult`` also carries a
        ``.units`` dict (str -> ``pint.Unit``) giving the unit of each key.

        Parameters
        ----------
        beam_energy : float
            Beam kinetic energy, in ``energy_unit``.
        angle_name : str
            Independent variable name, can be one of ``"theta3_lab"``, ``"theta4_lab"``,
            ``"theta_cm"``, ``"cos_theta_cm"``.
        angle_value : float
            Value to evaluate at. For angle keys (``theta*``), interpreted in
            ``angle_unit`` (default degrees). For ``"cos_theta_cm"``, treated
            as a dimensionless cosine — ``angle_unit`` is ignored.
        angle_unit : AngleUnit, optional
            Unit of ``angle_value`` for angle keys (default degrees).
        energy_unit : EnergyUnit, optional
            Unit of ``beam_energy`` (default MeV). Also governs the unit of every
            energy- and momentum-valued output and of ``duplicate_tol``, exactly
            as ``angle_unit`` governs every angle-valued output.
        duplicate_tol : float, optional
            Tolerance for merging near-duplicate solutions (default 1e-6), in
            ``energy_unit``. Solutions within this ``energy3_lab`` difference are
            treated as the same physical solution.

        Returns
        -------
        KinematicsResult[list[float]]
            Full dict of all kinematic variables, each a list of solutions sorted
            descending by ``energy3_lab``. Angle keys are in ``angle_unit``, energy
            and momentum keys in ``energy_unit`` (both default degrees/MeV) — see
            ``Reaction.output_units()`` for the full mapping, or just read
            ``.units[key]`` off the returned result directly.

        Raises
        ------
        ValueError
            If ``beam_energy`` or ``angle_value`` is not finite, if the reaction is
            kinematically forbidden at this energy, or if ``angle_value`` is outside
            the physical range.

        Examples
        --------
        >>> rxn = Reaction("p", "3H", "n", "3He")
        >>> rxn.kinematics_at_beam_energy_and_angle(1.2, "theta3_lab", 30)
        {'theta3_lab': [...], 'energy3_lab': [...], ...}
        """
        angle_unit = AngleUnit.from_any(angle_unit)
        energy_unit = EnergyUnit.from_any(energy_unit)
        output = self._kinematics_at_beam_energy_and_angle_raw(
            beam_energy,
            angle_name,
            angle_value,
            angle_unit=angle_unit,
            energy_unit=energy_unit,
            duplicate_tol=duplicate_tol,
        )
        return KinematicsResult(output, _result_units(output.keys(), angle_unit, energy_unit))

    def _kinematics_at_beam_energy_and_angle_raw(
        self,
        beam_energy: float,
        angle_name: str,
        angle_value: float,
        *,
        angle_unit: AngleUnit,
        energy_unit: EnergyUnit,
        duplicate_tol: float,
    ) -> dict[str, list[float]]:
        """Plain-float core of kinematics_at_beam_energy_and_angle, reused internally
        by kinematics_curve_at_angle."""
        # Keep the value/unit as the caller wrote it for error messages — the
        # internal, canonical-radians angle_value below is not what they typed.
        is_angle = angle_name.startswith("theta")
        angle_value_display = f"{angle_value} {angle_unit.name}" if is_angle else f"{angle_value}"
        if is_angle:
            angle_value = angle_value * angle_unit.value

        if not math.isfinite(angle_value):
            raise ValueError(f"{angle_name}={angle_value_display} is not a finite number")

        duplicate_tol_mev = duplicate_tol * energy_unit.value

        ek_mev = _parse_energy(beam_energy, energy_unit)
        self._bind(ek_mev)
        if self._nogo:
            raise ValueError(f"Reaction kinematically forbidden at beam_energy={ek_mev} MeV")

        if self._table is None:
            self._build_table()
        assert self._table is not None

        keys = list(self._table.keys())
        xs = self._table[angle_name]

        solutions = []

        exact_idx = np.where(np.isclose(xs, angle_value, atol=1e-12))[0]
        if len(exact_idx) > 0:
            for i in exact_idx:
                solutions.append({k: self._table[k][i] for k in keys})
        else:
            found = False
            for i in range(len(xs) - 1):
                x0, x1 = xs[i], xs[i + 1]
                if (x0 - angle_value) * (x1 - angle_value) <= 0 and x0 != x1:
                    found = True
                    t = (angle_value - x0) / (x1 - x0)
                    solutions.append(
                        {
                            k: self._table[k][i] + t * (self._table[k][i + 1] - self._table[k][i])
                            for k in keys
                        }
                    )
            if not found:
                raise ValueError(f"{angle_name}={angle_value_display} outside physical range")

        unique: list = []
        for sol in solutions:
            if not any(
                abs(sol["energy3_lab"] - u["energy3_lab"]) < duplicate_tol_mev for u in unique
            ):
                unique.append(sol)
        unique.sort(key=lambda s: s["energy3_lab"], reverse=True)

        output = {k: [s[k] for s in unique] for k in keys}
        for k in keys:
            if k in _ANGLE_KEYS:
                output[k] = [v / angle_unit.value for v in output[k]]
            elif k in _ENERGY_KEYS:
                output[k] = [v / energy_unit.value for v in output[k]]
        return output

    def kinematics_curve_at_angle(
        self,
        beam_energy_array: Iterable[float],
        theta3_lab: float,
        *,
        angle_unit: AngleUnit = AngleUnit.deg,
        energy_unit: EnergyUnit = EnergyUnit.MeV,
    ) -> list[KinematicsResult[npt.NDArray[np.float64]]]:
        """
        Compute ejectile kinematics at a fixed ejectile lab angle (``theta3_lab``)
        over a range of beam energies.

        Returns two branches (high- and low-energy) as a list of two
        ``KinematicsResult``s. Each contains plain arrays indexed by beam energy,
        with ``NaN`` where that branch does not exist, plus a ``.units`` dict
        (str -> ``pint.Unit``) computed from the same ``angle_unit``/``energy_unit``
        as the data. Branch 0 is always the higher-energy solution.

        Parameters
        ----------
        beam_energy_array : array-like
            Beam energies to sweep, in ``energy_unit``.
        theta3_lab : float
            Fixed ejectile lab angle (``theta3_lab`` in the output dict), in
            ``angle_unit``.
        angle_unit : AngleUnit, optional
            Unit of ``theta3_lab`` input and ``theta4_lab`` output (default degrees).
        energy_unit : EnergyUnit, optional
            Unit of ``beam_energy_array`` values (default MeV). Also governs the
            unit of every energy- and momentum-valued output (``beam_energy_lab``,
            ``energy3_lab``, ``energy4_lab``, ``momentum3_lab``, ``momentum4_lab``),
            exactly as ``angle_unit`` governs ``theta4_lab``.

        Returns
        -------
        list[KinematicsResult[np.ndarray]]
            List of two results, each with keys ``"beam_energy_lab"``, ``"energy3_lab"``,
            ``"energy4_lab"``, ``"theta4_lab"``, ``"velocity3_lab"``, ``"velocity4_lab"``,
            ``"momentum3_lab"``, ``"momentum4_lab"``, ``"jacobian3_lab"``, ``"jacobian4_lab"``.
            Angle keys are in ``angle_unit``, energy and momentum keys in
            ``energy_unit`` (both default degrees/MeV) — see ``Reaction.output_units()``
            for the full mapping, or just read ``.units[key]`` off either result.

        Examples
        --------
        >>> rxn = Reaction("p", "3H", "n", "3He")
        >>> branches = rxn.kinematics_curve_at_angle(np.linspace(1.0, 5.0, 200), 30)
        >>> for b in branches:
        ...     plt.plot(b["beam_energy_lab"], b["energy3_lab"])
        """
        angle_unit = AngleUnit.from_any(angle_unit)
        energy_unit = EnergyUnit.from_any(energy_unit)
        theta_rad = theta3_lab * angle_unit.value

        keys = [
            "energy3_lab",
            "energy4_lab",
            "theta4_lab",
            "velocity3_lab",
            "velocity4_lab",
            "momentum3_lab",
            "momentum4_lab",
            "jacobian3_lab",
            "jacobian4_lab",
        ]
        branches = [
            {"beam_energy_lab": [], **{k: [] for k in keys}},
            {"beam_energy_lab": [], **{k: [] for k in keys}},
        ]

        for ek in beam_energy_array:
            ek_mev = _parse_energy(ek, energy_unit)
            try:
                row = self._kinematics_at_beam_energy_and_angle_raw(
                    ek_mev,
                    "theta3_lab",
                    theta_rad,
                    angle_unit=AngleUnit.rad,
                    energy_unit=EnergyUnit.MeV,
                    duplicate_tol=1e-6,
                )
            except ValueError:
                solutions = []
            else:
                n = len(row["energy3_lab"])
                solutions = [{k: row[k][i] for k in keys} for i in range(n)]

            for i, branch in enumerate(branches):
                branch["beam_energy_lab"].append(ek_mev)
                sol = solutions[i] if i < len(solutions) else None
                for k in keys:
                    branch[k].append(sol[k] if sol is not None else float("nan"))

        result = []
        for branch in branches:
            converted = {}
            for k, v in branch.items():
                arr = np.array(v)
                if k in _ANGLE_KEYS:
                    arr = arr / angle_unit.value
                elif k in _ENERGY_KEYS:
                    arr = arr / energy_unit.value
                converted[k] = arr
            result.append(
                KinematicsResult(
                    converted, _result_units(converted.keys(), angle_unit, energy_unit)
                )
            )
        return result

    @staticmethod
    def output_units(
        *,
        angle_unit: AngleUnit | str = AngleUnit.deg,
        energy_unit: EnergyUnit | str = EnergyUnit.MeV,
    ) -> dict[str, pint.Unit]:
        """
        The ``pint.Unit`` for every key that can appear in the dicts returned by
        ``kinematics_table_at_beam_energy``, ``kinematics_at_beam_energy_and_angle``,
        and ``kinematics_curve_at_angle`` — the same mapping each of those methods
        attaches as ``result.units``, available here without running a computation
        first.

        Pass the same ``angle_unit``/``energy_unit`` you'd pass to one of those
        methods to get a matching ``{key: pint.Unit}`` map — not every key appears
        in every method's output, so look up only the keys actually present in the
        result you're inspecting. ``velocity3_lab``/``velocity4_lab`` (fraction of
        the speed of light) and ``jacobian3_lab``/``jacobian4_lab`` (ratio of solid
        angles, dOmega_lab/dOmega_cm) are both genuinely dimensionless, not just
        unlabeled.

        Examples
        --------
        >>> str(Reaction.output_units()["theta3_lab"])
        'degree'
        >>> f"{Reaction.output_units(energy_unit='keV')['momentum3_lab']:~}"
        'keV / c'
        """
        angle_unit = AngleUnit.from_any(angle_unit)
        energy_unit = EnergyUnit.from_any(energy_unit)
        # Built from the same _result_unit used to actually tag KinematicsResult.units,
        # so this preview can't drift out of sync with real results.
        return _result_units(_ALL_RESULT_KEYS, angle_unit, energy_unit)
