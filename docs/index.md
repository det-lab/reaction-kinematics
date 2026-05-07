# Welcome to Reaction Kinematics 
This is a Python library for calculating relativistic two-body nuclear reaction kinematics.

This package is designed for students and researchers working in nuclear and particle physics who need fast, reliable kinematic calculations for reactions of the form:

```
projectile + target → ejectile + recoil
```

This code can do:

 * Relativistic two-body kinematics
 * Automatic unit handling
 * Center-of-mass and lab-frame quantities
 * Energy, angle, momentum, and velocity calculations
 * Support for multi-valued kinematic solutions
 * Simple plotting and data export
<br>

## About
The main interface is the `Reaction` class.

```python
from reaction_kinematics import Reaction
```

Create a reaction by specifying the four particle masses. The beam energy is passed separately to each calculation method.

### Reaction 
We define a reaction with these variables

* m1 (Projectile)
* m2 (Target)
* m3 (Ejectile)
* m4 (Recoil)

All masses are converted internally to MeV. Kinetic energy is assumed in MeV unless specified otherwise.

For example if you wanted to do p(3H, n)3He:
```python
from reaction_kinematics import Reaction

rxn = Reaction("p", "3H", "n", "3He")
```
<br> 

### Units 

* Masses are internally stored in MeV/c²
* Energies are in MeV by default
* Velocities are given as fractions of c
* Angles are in radians
<br>

### Compute Arrays


To generate arrays of kinematic quantities over all center-of-mass angles, use `kinematics_table_at_beam_energy()`.

```python
data = rxn.kinematics_table_at_beam_energy(4.0)
```

This will return a dictionary containing the following:

* `coscm`   : cos(θ_CM)
* `theta_cm`: CM angle (rad)
* `theta3`  : Ejectile lab angle (rad)
* `theta4`  : Recoil lab angle (rad)
* `e3`      : Ejectile energy (MeV)
* `e4`      : Recoil energy (MeV)
* `v3`      : Ejectile velocity (c)
* `v4`      : Recoil velocity (c)
<br>

###### See [Plotting Example](plotting.md) for How to Plot


