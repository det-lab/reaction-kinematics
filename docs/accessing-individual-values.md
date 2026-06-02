## Accessing Individual Values

To evaluate kinematic quantities at a specific beam energy and kinematic variable value, use `kinematics_at_beam_energy_and_angle()`.

This method automatically handles multi-valued solutions and always returns lists.

### Syntax

```python
from reaction_kinematics import Reaction

# or equivalently: Reaction("3H(p,n)3He")
rxn = Reaction("p", "3H", "n", "3He")

r = rxn.kinematics_at_beam_energy_and_angle(1.2, "theta3_lab", 122)
print(r)

```
This will return:
```python
{'cos_theta_cm': [-0.9813047694699246], 
'theta_cm': [168.9], 
'theta3_lab': [122.0], 
'theta4_lab': [3.098], 
'energy3_lab': [0.005237930117122223], 
'energy4_lab': [0.43100601988343396], 
'velocity3_lab': [0.003338681989975168], 
'velocity4_lab': [0.017517730726598162], 
'momentum3_lab': [3.136927647029636], 
'momentum4_lab': [49.204196838704895]}
```
Parameters:

* `beam_energy`  : Beam kinetic energy (MeV by default)
* `angle_name`   : Independent variable (e.g. `"theta3_lab"`, `"theta4_lab"`, `"theta_cm"`, `"cos_theta_cm"`)
* `angle_value`  : Value at which to evaluate (degrees for angles by default)
* `angle_unit`   : Unit for `angle_value` — `"deg"` (default), `"rad"`, `"mrad"`
* `energy_unit`  : Unit for `beam_energy` — `"keV"`, `"MeV"` (default), `"GeV"`, `"TeV"`

For example, to evaluate at a beam energy given in keV:
```python
r = rxn.kinematics_at_beam_energy_and_angle(1200, "theta3_lab", 122, energy_unit="keV")
```
