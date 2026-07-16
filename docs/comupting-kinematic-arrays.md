## Computing Kinematic Arrays

To generate arrays of kinematic quantities over all center-of-mass angles, use `kinematics_table_at_beam_energy()`.

```python
from reaction_kinematics import Reaction

# or equivalently: Reaction("3H(p,n)3He")
rxn = Reaction("p", "3H", "n", "3He")

data = rxn.kinematics_table_at_beam_energy(1.2)

print(data["theta4_lab"])
print(data["energy3_lab"])
```

The `energy_unit` parameter accepts `"keV"`, `"MeV"` (default), `"GeV"`, or `"TeV"`.
It governs `energy3_lab`, `energy4_lab`, `momentum3_lab`, and `momentum4_lab` in
the returned result too, not just the input beam energy:
```python
data = rxn.kinematics_table_at_beam_energy(1200, energy_unit="keV")  # equivalent to 1.2 MeV
# energy3_lab, energy4_lab, momentum3_lab, momentum4_lab are all in keV here
```

This will return a `KinematicsResult` (behaves like a dict: `data["energy3_lab"]`)
containing the following:

* `cos_theta_cm`   : cos(θ_CM)
* `theta_cm`       : CM angle (deg)
* `theta3_lab`     : Ejectile lab angle (deg)
* `theta4_lab`     : Recoil lab angle (deg)
* `energy3_lab`    : Ejectile energy (MeV)
* `energy4_lab`    : Recoil energy (MeV)
* `velocity3_lab`  : Ejectile velocity (unitless, reported as a fraction of c)
* `velocity4_lab`  : Recoil velocity (unitless, reported as a fraction of c)
* `momentum3_lab`  : Ejectile momentum (MeV/c)
* `momentum4_lab`  : Recoil momentum (MeV/c)
* `jacobian3_lab`  : dΩ₃(lab)/dΩ(cm), converts a lab-frame differential cross section to the cm frame
* `jacobian4_lab`  : dΩ₄(lab)/dΩ(cm), converts a lab-frame differential cross section to the cm frame

Each key's unit: `data.units["energy3_lab"]`, or `Reaction.output_units(energy_unit="keV")` for a preview without running a computation.


