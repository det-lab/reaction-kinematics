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

For example, for the reaction ³H(p,n)³He:
```python
from reaction_kinematics import Reaction

# or equivalently: Reaction("3H(p,n)3He")
rxn = Reaction("p", "3H", "n", "3He")
```
<br> 

### Units 

* Masses are internally stored in MeV/c²
* Energies are in MeV by default — supported units: `keV`, `MeV`, `GeV`, `TeV`
* Velocities are given as fractions of c
* Angles are in degrees by default — supported units: `deg`, `rad`, `mrad`
<br>

### Compute Arrays


To generate arrays of kinematic quantities over all center-of-mass angles, use `kinematics_table_at_beam_energy()`.

```python
data = rxn.kinematics_table_at_beam_energy(4.0)
```

This will return a dictionary containing the following:

* `cos_theta_cm`  : cos(θ_CM)
* `theta_cm`      : CM angle (deg)
* `theta3_lab`    : Ejectile lab angle (deg)
* `theta4_lab`    : Recoil lab angle (deg)
* `energy3_lab`   : Ejectile energy (MeV)
* `energy4_lab`   : Recoil energy (MeV)
* `velocity3_lab` : Ejectile velocity (c)
* `velocity4_lab` : Recoil velocity (c)
* `momentum3_lab` : Ejectile momentum (MeV/c)
* `momentum4_lab` : Recoil momentum (MeV/c)
<br>

###### See [Plotting Example](plotting.md) for How to Plot


