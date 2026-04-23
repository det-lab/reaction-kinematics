## Computing Kinematic Arrays

To generate arrays of kinematic quantities over all center-of-mass angles, use `compute_arrays()`.

```python
import numpy as np
from reaction_kinematics import TwoBody

rxn = TwoBody("p", "3H", "n", "3He", 1.2)
#Proton + Tritium Reaction

data = rxn.compute_arrays()

theta4 = np.array(data["theta4"])
e3 = np.array(data["e3"])

print(theta4)
print(e3)
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



