# reaction-kinematics
This is a Python library for calculating relativistic two-body nuclear reaction kinematics.

This package is designed for students and researchers working in nuclear and particle physics who need reliable relativistic kinematic calculations for reactions of the form:

```
projectile + target → ejectile + recoil
```

---

## Features 
This code can do:
* Relativistic two-body kinematics
* Center-of-mass and lab-frame quantities
* Energy, angle, momentum, and velocity calculations
* Support for multi-valued kinematic solutions

---

## Installation


```
pip install reaction-kinematics
```


---

## Documentation 

https://det-lab.github.io/reaction-kinematics/


### Syntax

```python
# or equivalently:
# rxn = Reaction("3H(p,n)3He")
rxn = Reaction("p", "3H", "n", "3He")

rxn.kinematics_at_beam_energy_and_angle(beam_energy, angle_name, angle_value)
```

Parameters:

* `beam_energy`  : Beam kinetic energy (MeV by default)
* `angle_name`   : Independent variable (e.g. `"theta3_lab"`, `"theta4_lab"`, `"theta_cm"`, `"cos_theta_cm"`)
* `angle_value`  : Value at which to evaluate (degrees for angles by default)
* `angle_unit`   : Unit for `angle_value` — `"deg"` (default), `"rad"`, `"mrad"`
* `energy_unit`  : Unit for `beam_energy` — `"keV"`, `"MeV"` (default), `"GeV"`, `"TeV"`

---

## Example

```python
vals = rxn.kinematics_at_beam_energy_and_angle(1.2, "theta4_lab", 10)
```

`vals` is a `KinematicsResult` — behaves like a dict (`vals["energy3_lab"]`), plus
a `.units` dict giving the `pint.Unit` of each key (`vals.units["energy3_lab"]`):

```
{
  'cos_theta_cm': [...],
  'theta_cm': [...],
  'theta3_lab': [...],
  'theta4_lab': [...],
  'energy3_lab': [0.3448, 0.0364],
  'energy4_lab': [...],
  'velocity3_lab': [0.027, 0.009],
  'velocity4_lab': [...],
  'momentum3_lab': [25.5, 8.3],
  'momentum4_lab': [...],
  'jacobian3_lab': [...],
  'jacobian4_lab': [...]
}
```

Multiple values in each list indicate multiple physical solutions.
`jacobian3_lab`/`jacobian4_lab` are dΩ(lab)/dΩ(cm) for the ejectile/recoil —
the factor that converts a lab-frame differential cross section to the cm frame.

---

### Using Explicit Mass Values

```python
rxn = Reaction(
    938.272,
    11177.928,
    938.272,
    11177.928,
    mass_unit="MeV"
)
```

---


## Numerical Notes

* Some kinematic variables are multi-valued.
* Near kinematic extrema, solution branches may merge numerically.
* The library automatically removes duplicate solutions within tolerance of 1e-6.

---

## License

GPL-2.0 license

---

## Contact

For questions, issues, or contributions, please open an issue on GitHub.


* * *

## Project Docs

For how to install uv and Python, see [installation.md](installation.md).

For development workflows, see [development.md](development.md).

For instructions on publishing to PyPI, see [publishing.md](publishing.md).

* * *

*This project was built from
[simple-modern-uv](https://github.com/jlevy/simple-modern-uv).*

