"""
Relativistic two-body reaction kinematics
"""

import math

import numpy as np

from reaction_kinematics.inputs import MassInput
from reaction_kinematics.units import AngleUnit, EnergyUnit


def _parse_mass(m, unit=None):
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


def q_value(m1, m2, m3, m4, *, mass_unit=None) -> float:
    """
    Compute the Q-value for a reaction: projectile + target → ejectile + recoil

    Q = m1 + m2 - m3 - m4  (in MeV)

    Accepts the same mass input types as TwoBody: string symbols, MassInput objects,
    or numeric values (which require mass_unit).

    Parameters
    ----------
    m1, m2, m3, m4 : str, MassInput, or float
        Masses of projectile, target, ejectile, and recoil.
    mass_unit : str or EnergyUnit, optional
        Unit for numeric masses (e.g. "MeV", "keV", "amu").

    Returns
    -------
    float
        Q-value in MeV.
    """
    return (
        _parse_mass(m1, mass_unit)
        + _parse_mass(m2, mass_unit)
        - _parse_mass(m3, mass_unit)
        - _parse_mass(m4, mass_unit)
    )


def kinematic_curve(m1, m2, m3, m4, theta3, ek_array, *, mass_unit=None, angle_unit=AngleUnit.rad):
    """
    Compute ejectile kinematics at a fixed lab angle over a range of beam energies.

    Returns two branches (high-energy and low-energy) as a list of two dicts. Each
    dict contains arrays indexed by beam energy, with NaN where that branch does not
    exist. Branch 0 is always the higher-energy solution.

    Parameters
    ----------
    m1, m2, m3, m4 : str, MassInput, or float
        Masses of projectile, target, ejectile, and recoil.
    theta3 : float
        Fixed lab angle of the ejectile.
    ek_array : array-like
        Beam energies to evaluate (in MeV unless energy_unit is specified via mass_unit).
    mass_unit : str or EnergyUnit, optional
        Unit for numeric masses.
    angle_unit : AngleUnit or str, optional
        Unit of theta3 (default: radians).

    Returns
    -------
    list of two dicts, each with keys:
        - "ek"     : beam energy (MeV)
        - "e3"     : ejectile kinetic energy (MeV)
        - "e4"     : recoil kinetic energy (MeV)
        - "theta4" : recoil lab angle (radians)
        - "v3"     : ejectile velocity (fraction of c)
        - "v4"     : recoil velocity (fraction of c)

    Example
    -------
    >>> import numpy as np
    >>> branches = kinematic_curve("p", "3H", "n", "3He", np.deg2rad(30), np.linspace(1.0, 5.0, 200))
    >>> for b in branches:
    ...     plt.plot(b["ek"], b["e3"])
    """
    if isinstance(angle_unit, str):
        angle_unit = AngleUnit[angle_unit]

    theta3_rad = theta3 * angle_unit.value

    keys = ["e3", "e4", "theta4", "v3", "v4"]
    branches = [
        {"ek": [], **{k: [] for k in keys}},
        {"ek": [], **{k: [] for k in keys}},
    ]

    for ek in ek_array:
        rxn = TwoBody(m1, m2, m3, m4, ek, mass_unit=mass_unit)

        try:
            row = rxn.at_value("theta3", theta3_rad, y_names=keys)
        except ValueError:
            solutions = []
        else:
            n = len(row["e3"])
            solutions = [{k: row[k][i] for k in keys} for i in range(n)]

        for i, branch in enumerate(branches):
            branch["ek"].append(ek)
            sol = solutions[i] if i < len(solutions) else None
            for k in keys:
                branch[k].append(sol[k] if sol is not None else float("nan"))

    return [{k: np.array(v) for k, v in branch.items()} for branch in branches]


class TwoBody:
    """
    This class calculates the lab-frame and center-of-mass kinematics for a reaction:
        projectile + target → ejectile + recoil

    It automatically computes:
        - Center-of-mass energies and momenta
        - Lab-frame energy extrema
        - Lab-frame angles and velocities

    All masses are converted internally to MeV. Kinetic energy is assumed in MeV unless specified otherwise.
        Angles are returned in radians by default. Use the `AngleUnit` enum to convert if needed.

    where:
        m1 : Projectile (MassInput, string symbol ["alpha"], or numeric mass [ 3727.3794118])
        m2 : Target (MassInput, string symbol["12C], or numeric mass [11177.928])
        m3 : Ejectile (MassInput, string symbol ["p"], or numeric mass [938.272])
        m4 : Recoil (MassInput, string symbol ["12C"], or numeric mass[11177.928])

    Parameters
    ----------
    m1, m2, m3, m4 : str, MassInput, or float
        Masses of the particles. Strings like "p", "12C" will be converted
        using the internal mass table. Floats must be accompanied by `mass_unit`.
    ek : float or EnergyValue
        Kinetic energy of the projectile in the lab frame. If a float is given,
        `energy_unit` specifies the units (default: MeV).
    mass_unit : str or EnergyUnit, optional
        Unit of numeric masses if floats are provided (e.g., "MeV", "keV", "amu").
    energy_unit : EnergyUnit or str, optional
        Unit of kinetic energy (default: MeV).

    Attributes
    ----------
    s : float
        Mandelstam s of the reaction.
    pcm : float
        Initial center-of-mass momentum.
    pcmp : float
        Final-state center-of-mass momentum.
    e03, e04 : float
        Center-of-mass energies of particles 3 and 4.
    emax3, emin3, emax4, emin4 : float
        Lab-frame extrema of particle energies.
    theta3max, theta4max : float
        Maximum lab-frame angles for ejectile and recoil (radians).
    cmcos3max, cmcos4max : float
        Maximum CM frame cos(theta) for particles 3 and 4.

    User-facing methods
    ------------------
    compute_arrays()
        Returns a dictionary of arrays for cos(theta_CM), theta3, theta4, e3, e4, v3, v4.
    get_points()
        Returns a list of dictionaries with kinematic values at each CM angle.
    at_value(quantity, theta_cm)
        Returns interpolated values of a chosen quantity at a given CM angle.

    Returned dictionary from compute_arrays():
        - "coscm"   : cos(theta) in the CM frame
        - "theta_cm": theta in the CM frame (radians)
        - "theta3"  : lab-frame angle of ejectile (radians)
        - "theta4"  : lab-frame angle of recoil (radians)
        - "e3"      : lab-frame kinetic energy of ejectile (MeV)
        - "e4"      : lab-frame kinetic energy of recoil (MeV)
        - "v3"      : lab-frame velocity of ejectile (fraction of c)
        - "v4"      : lab-frame velocity of recoil (fraction of c)

    Example
    -------
    # Using MassInput objects
    >>> from reaction_kinematics.inputs import MassInput
    >>> p = MassInput("p")       # proton
    >>> C12 = MassInput("12C")   # carbon-12
    >>> ek = 5.0                 # MeV
    >>> rxn = TwoBody(p, C12, p, C12, ek)
    >>> data = rxn.compute_arrays()
    >>> print(data["theta3"][:5])

    # Using string symbols (auto-converted using mass table)
    >>> rxn2 = TwoBody("p", "12C", "p", "12C", ek)
    >>> data2 = rxn2.compute_arrays()
    >>> print(data2["e3"][:5])

    # Using numeric masses with a specified unit
    >>> # masses in MeV
    >>> m_proj = 938.272    # proton
    >>> m_targ = 11177.928  # carbon-12
    >>> m_ej = 938.272
    >>> m_recoil = 11177.928
    >>> rxn3 = TwoBody(m_proj, m_targ, m_ej, m_recoil, ek, mass_unit="MeV")
    >>> data3 = rxn3.compute_arrays()
    >>> print(data3["v3"][:5])

    # Quick start:
    >>>  from reaction_kinematics.reaction_kinematics import TwoBody
    >>> rxn = TwoBody("p", "12C", "p", "12C", 5.0)
    >>> data = rxn.compute_arrays()
    """

    def __init__(self, m1, m2, m3, m4, ek, *, mass_unit=None, energy_unit=EnergyUnit.MeV):
        self.m1 = _parse_mass(m1, mass_unit)
        self.m2 = _parse_mass(m2, mass_unit)
        self.m3 = _parse_mass(m3, mass_unit)
        self.m4 = _parse_mass(m4, mass_unit)

        if isinstance(ek, (int, float)):
            if isinstance(energy_unit, str):
                energy_unit = EnergyUnit[energy_unit]
            self.ek = ek * energy_unit.value
        else:
            self.ek = ek

        # defaults
        self.ncoscm = 500
        self.nogo = False

        self.cmcos3max = 2.0
        self.cmcos4max = 2.0
        self.e3atmaxang = -1.0
        self.e4atmaxang = -1.0
        self.theta3max = None
        self.theta4max = None

        # kinematic quantities (always defined)
        self.s = None
        self.pcm = None
        self.pcmp = None
        self.cmrap = None
        self.thesinh = None
        self.thecosh = None
        self.e03 = None
        self.e04 = None
        self.emax3 = None
        self.emin3 = None
        self.emax4 = None
        self.emin4 = None

        self._compute()

    def _compute(self):
        # Mandelstam s
        self.s = (self.m1 + self.m2) ** 2 + 2.0 * self.m2 * self.ek
        if self.s <= 0.0:
            self.nogo = True
            return

        # initial CM momentum
        pcm2 = (self.s - self.m1**2 - self.m2**2) ** 2 - 4.0 * self.m1**2 * self.m2**2
        if pcm2 < 0:
            self.nogo = True
            return

        self.pcm = math.sqrt(pcm2 / (4.0 * self.s))

        # CM rapidity
        acmratio = (math.sqrt(self.m2**2 + self.pcm**2) + self.pcm) / self.m2
        self.cmrap = math.log(acmratio)
        self.thesinh = math.sinh(self.cmrap)
        self.thecosh = math.cosh(self.cmrap)

        # final-state CM momentum
        pcmp2 = (self.s - self.m3**2 - self.m4**2) ** 2 - 4.0 * self.m3**2 * self.m4**2
        if pcmp2 < 0:
            self.nogo = True
            return

        self.pcmp = math.sqrt(pcmp2 / (4.0 * self.s))

        # CM energies
        self.e03 = math.sqrt(self.pcmp**2 + self.m3**2)
        self.e04 = math.sqrt(self.pcmp**2 + self.m4**2)

        # lab-frame extrema
        self.emax3 = self.e03 * self.thecosh + self.pcmp * self.thesinh - self.m3
        self.emin3 = self.e03 * self.thecosh - self.pcmp * self.thesinh - self.m3
        self.emax4 = self.e04 * self.thecosh + self.pcmp * self.thesinh - self.m4
        self.emin4 = self.e04 * self.thecosh - self.pcmp * self.thesinh - self.m4

        # max ejectile angle
        thetatest: float | None = None
        if self.m3 > 0.0:
            thetatest = self.pcmp / (self.m3 * self.thesinh)

            if thetatest is not None and thetatest < 1.0:
                self.theta3max = math.asin(thetatest)

                patmax = (self.e03 * math.cos(self.theta3max) * self.thesinh) / (
                    1.0 + thetatest**2 * self.thesinh**2
                )

                eatmax = math.sqrt(patmax**2 + self.m3**2)
                self.e3atmaxang = eatmax - self.m3

                self.cmcos3max = (eatmax - self.e03 * self.thecosh) / (self.pcmp * self.thesinh)

        # elastic symmetry case
        if (self.m1 + self.m2) == (self.m3 + self.m4):
            if thetatest is None:
                raise RuntimeError("thetatest was not computed")
            if abs(thetatest - 1.0) < 1e-3:
                self.theta3max = math.pi / 2.0
                self.cmcos3max = -1.0
                self.e3atmaxang = (
                    self.e03 * self.thecosh + self.cmcos3max * self.pcmp * self.thesinh - self.m3
                )

        # max recoil angle
        if self.m4 > 0.0:
            thetatest = self.pcmp / (self.m4 * self.thesinh)
            if thetatest is not None and thetatest < 1.0:
                self.theta4max = math.asin(thetatest)
                patmax = (self.e04 * math.cos(self.theta4max) * self.thesinh) / (
                    1.0 + thetatest**2 * self.thesinh**2
                )
                eatmax = math.sqrt(patmax**2 + self.m4**2)
                self.e4atmaxang = eatmax - self.m4
                self.cmcos4max = (eatmax - self.e04 * self.thecosh) / (self.pcmp * self.thesinh)

        if (self.m1 + self.m2) == (self.m3 + self.m4):
            if thetatest is None:
                raise RuntimeError("thetatest was not computed")
            if abs(thetatest - 1.0) < 1e-3:
                self.theta4max = math.pi / 2.0
                self.cmcos4max = 1.0
                self.e4atmaxang = (
                    self.e04 * self.thecosh - self.cmcos4max * self.pcmp * self.thesinh - self.m4
                )

    @property
    def q_value(self) -> float:
        """Q-value of the reaction in MeV: Q = m1 + m2 - m3 - m4"""
        return self.m1 + self.m2 - self.m3 - self.m4

    def _kinematics_at_coscm(self, coscm):
        # Kinematic computations require these to have been initialized by _compute()
        # For example if reaction is below threshold these values will not be computed
        if (
            self.pcmp is None
            or self.thecosh is None
            or self.e03 is None
            or self.thesinh is None
            or self.e04 is None
        ):
            raise ValueError("Kinematic quantities not computed")
        sincm = math.sqrt(max(0.0, 1.0 - coscm**2))

        ppar3 = self.pcmp * self.thecosh * coscm + self.e03 * self.thesinh
        pperp3 = self.pcmp * sincm
        ptot3 = math.hypot(ppar3, pperp3)

        ppar4 = -self.pcmp * self.thecosh * coscm + self.e04 * self.thesinh
        pperp4 = self.pcmp * sincm
        ptot4 = math.hypot(ppar4, pperp4)

        theta3 = math.acos(ppar3 / ptot3) if ptot3 > 0 else 0.0
        theta4 = math.acos(ppar4 / ptot4) if ptot4 > 0 else 0.0

        e3 = self.e03 * self.thecosh + coscm * self.pcmp * self.thesinh - self.m3
        e4 = self.e04 * self.thecosh - coscm * self.pcmp * self.thesinh - self.m4

        v3 = ptot3 / (e3 + self.m3)
        v4 = ptot4 / (e4 + self.m4)

        return {
            "coscm": coscm,
            "theta_cm": math.acos(coscm),
            "theta3": theta3,
            "theta4": theta4,
            "e3": e3,
            "e4": e4,
            "v3": v3,
            "v4": v4,
            "p3": ptot3,
            "p4": ptot4,
        }

    def _solve_at_theta_cm(self, theta_cm):
        coscm = math.cos(theta_cm)
        return self._kinematics_at_coscm(coscm)

    def _build_table(self):
        table = {}
        keys = ["coscm", "theta_cm", "theta3", "theta4", "e3", "e4", "v3", "v4", "p3", "p4"]
        for k in keys:
            table[k] = []

        for i in range(-self.ncoscm, self.ncoscm + 1):
            coscm = i / self.ncoscm
            row = self._kinematics_at_coscm(coscm)
            for k in keys:
                table[k].append(row[k])

        self._table = table

    def at_value(self, x_name, x, *, y_names=None, duplicate_tol=1e-6):
        """
        Interpolate kinematic quantities.

        Always returns lists to handle multi-valued cases.

        Parameters
        ----------
        x_name : str
            Independent variable (e.g. 'theta3', 'theta_cm', 'coscm')
        x : float
            Value to evaluate at
        y_names : list[str] or None
            Dependent variables. If None, returns all.
        duplicate_tol: float
            A point where it will stop two roots and merge into one
        """

        if not hasattr(self, "_table"):
            self._build_table()

        xs = self._table[x_name]

        if isinstance(y_names, str):
            y_names = [y_names]

        if y_names is None:
            y_names = list(self._table.keys())

        # Collect solutions as tuples to keep keys in sync across deduplication/sorting.
        solutions = []

        exact_idx = np.where(np.isclose(xs, x, atol=1e-12))[0]
        if len(exact_idx) > 0:
            for i in exact_idx:
                solutions.append({k: self._table[k][i] for k in y_names})
        else:
            # Interpolation
            found = False
            for i in range(len(xs) - 1):
                x0, x1 = xs[i], xs[i + 1]
                if (x0 - x) * (x1 - x) <= 0 and x0 != x1:
                    found = True
                    t = (x - x0) / (x1 - x0)
                    sol = {}
                    for k in y_names:
                        y0 = self._table[k][i]
                        y1 = self._table[k][i + 1]
                        sol[k] = y0 + t * (y1 - y0)
                    solutions.append(sol)

            if not found:
                raise ValueError(f"{x_name}={x} outside physical range")

        # Deduplicate on e3 (or first available key), then sort descending by e3.
        e_key = "e3" if "e3" in y_names else y_names[0]
        unique_solutions = []
        for sol in solutions:
            if not any(abs(sol[e_key] - u[e_key]) < duplicate_tol for u in unique_solutions):
                unique_solutions.append(sol)
        unique_solutions.sort(key=lambda s: s[e_key], reverse=True)

        results = {k: [s[k] for s in unique_solutions] for k in y_names}
        return results

    def _unique(self, arr, tol=1e-6):
        out = []
        for v in arr:
            if not any(abs(v - u) < tol for u in out):
                out.append(v)
        return out

    def compute_arrays(self):
         """
         Compute full two-body kinematics over center-of-mass angles.

         This method evaluates the reaction kinematics across a grid of
         center-of-mass (CM) angles and returns the corresponding lab-frame
         quantities for both outgoing particles.

         The results are stored internally and also returned as a dictionary
         of NumPy arrays. These arrays form the basis for interpolation via
         `at_value()`.

         Returns
         -------
         dict of str -> numpy.ndarray
             Dictionary containing kinematic quantities evaluated over the
             CM angle grid. Keys include:

             - "coscm"     : cos(theta_cm)
             - "theta_cm"  : center-of-mass angle (radians)
             - "theta3"    : ejectile lab angle (radians)
             - "theta4"    : recoil lab angle (radians)
             - "e3"        : ejectile kinetic energy (MeV)
             - "e4"        : recoil kinetic energy (MeV)
             - "p3"        : ejectile momentum
             - "p4"        : recoil momentum
             - "v3"        : ejectile velocity (fraction of c)
             - "v4"        : recoil velocity (fraction of c)

         Notes
         -----
         - The angular grid is defined internally using `self.ncoscm`.
         - All angles are in radians.
         - Energies are in MeV.
         - If the reaction is kinematically forbidden (e.g., below threshold),
           a ValueError is raised.

         Examples
         --------
         >>> rxn = TwoBody("p", "3H", "n", "3He", 1.0)
         >>> data = rxn.compute_arrays()
         >>> data["theta3"]  # array of lab angles

         This output can be used for plotting or passed to `at_value()` for
         interpolation at specific angles. 
         """
        data = {k: [] for k in ["coscm", "theta_cm", "theta3", "theta4", "e3", "e4", "v3", "v4"]}

        for i in range(-self.ncoscm, self.ncoscm + 1):
            coscm = i / self.ncoscm
            kin = self._kinematics_at_coscm(coscm)

            for k in data:
                data[k].append(kin[k])

        return data
