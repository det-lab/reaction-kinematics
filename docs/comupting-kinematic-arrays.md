## Computing Kinematic Arrays

To generate arrays of kinematic quantities over all center-of-mass angles, use `compute_arrays()`.

```python
data = rxn.compute_arrays()
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

### Example

```python
theta4 = data["theta4"]
e3 = data["e3"]
```

---
