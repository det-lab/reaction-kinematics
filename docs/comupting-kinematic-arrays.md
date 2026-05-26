## Computing Kinematic Arrays

To generate arrays of kinematic quantities over all center-of-mass angles, use `kinematics_table_at_beam_energy()`.

```python
from reaction_kinematics import Reaction

rxn = Reaction("p", "3H", "n", "3He")
# Proton + Tritium Reaction

data = rxn.kinematics_table_at_beam_energy(1.2)

print(data["theta4_lab"])
print(data["energy3_lab"])
```

This will return a dictionary containing the following:

* `cos_theta_cm`   : cos(θ_CM)
* `theta_cm`       : CM angle (rad)
* `theta3_lab`     : Ejectile lab angle (rad)
* `theta4_lab`     : Recoil lab angle (rad)
* `energy3_lab`    : Ejectile energy (MeV)
* `energy4_lab`    : Recoil energy (MeV)
* `velocity3_lab`  : Ejectile velocity (c)
* `velocity4_lab`  : Recoil velocity (c)


