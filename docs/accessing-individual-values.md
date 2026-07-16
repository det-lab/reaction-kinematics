## Accessing Individual Values

To evaluate kinematic quantities at a specific beam energy and kinematic variable value, use `kinematics_at_beam_energy_and_angle()`.

This method automatically handles multi-valued solutions and always returns lists.

### Syntax

```python
from reaction_kinematics import Reaction

# or equivalently: Reaction("3H(p,n)3He")
rxn = Reaction("p", "3H", "n", "3He")

r = rxn.kinematics_at_beam_energy_and_angle(1.2, "theta3_lab", 122)
print(r["energy3_lab"], r["theta4_lab"])
```
This returns a `KinematicsResult`, which behaves like a dict with these keys:

* `cos_theta_cm`: `[-0.981259434287812]`
* `theta_cm`: `[168.89373250634032]`
* `theta3_lab`: `[121.99999999999999]`
* `theta4_lab`: `[3.1017567808349322]`
* `energy3_lab`: `[0.005245865820720845]`
* `energy4_lab`: `[0.43099808417983787]`
* `velocity3_lab`: `[0.003341198974529349]`
* `velocity4_lab`: `[0.017517569485073135]`
* `momentum3_lab`: `[3.1392925586297746]`
* `momentum4_lab`: `[49.2037438008774]`
* `jacobian3_lab`: `[13.26018747199885]` — dΩ(lab)/dΩ(cm) for the ejectile, converts a lab-frame differential cross section to the cm frame
* `jacobian4_lab`: `[-0.07813620376523883]` — same, for the recoil

Each key's unit: `r.units["energy3_lab"]`.

Parameters:

* `beam_energy`  : Beam kinetic energy (MeV by default)
* `angle_name`   : Independent variable (e.g. `"theta3_lab"`, `"theta4_lab"`, `"theta_cm"`, `"cos_theta_cm"`)
* `angle_value`  : Value at which to evaluate (degrees for angles by default)
* `angle_unit`   : Unit for `angle_value` — `"deg"` (default), `"rad"`, `"mrad"`
* `energy_unit`  : Unit for `beam_energy` — `"keV"`, `"MeV"` (default), `"GeV"`, `"TeV"`

`energy_unit` also governs `energy3_lab`, `energy4_lab`, `momentum3_lab`, and
`momentum4_lab` in the result, not just the input beam energy. For example, to
evaluate at a beam energy given in keV:
```python
r = rxn.kinematics_at_beam_energy_and_angle(1200, "theta3_lab", 122, energy_unit="keV")
# energy3_lab, energy4_lab, momentum3_lab, momentum4_lab are all in keV here
```
