# reaction-kinematics
This is a Python library for calculating relativistic two-body nuclear reaction kinematics.

This package is designed for students and researchers working in nuclear and particle physics who need fast, reliable kinematic calculations for reactions of the form:

```
projectile + target → ejectile + recoil
```

---

## Features 
This code can do:
* Relativistic two-body kinematics
* Automatic unit handling
* Center-of-mass and lab-frame quantities
* Energy, angle, momentum, and velocity calculations
* Support for multi-valued kinematic solutions
* Simple plotting and data export

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
rxn.kinematics_at_beam_energy_and_variable(beam_energy, var_name, var_value, return_variables=None)
```

Parameters:

* `beam_energy` : Beam kinetic energy in MeV
* `var_name`    : Independent variable (e.g. `"theta4"`, `"theta_cm"`, `"coscm"`)
* `var_value`   : Value at which to evaluate (radians for angles)
* `return_variables`     : Dependent variables (string or list, `None` returns all)

---

## Example: Single Quantity

```python
import math

angle = 10 * math.pi / 180

vals = rxn.kinematics_at_beam_energy_and_variable(1.2, "theta4", angle, return_variables="e3")
print(vals)
```

Output:

```
{'e3': [0.3447, 0.0364]}
```

Multiple values indicate multiple physical solutions.

---

## Example: Multiple Quantities

```python
vals = rxn.kinematics_at_beam_energy_and_variable(
    1.2,
    "theta4",
    angle,
    return_variables=["e3", "v3", "p3"]
)

print(vals)
```

Example output:

```
{
  'e3': [0.3447, 0.0364],
  'v3': [0.025, 0.009],
  'p3': [23.7, 8.2]
}
```

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
* The library automatically removes duplicate solutions within tolerance of 1e**-6.

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

