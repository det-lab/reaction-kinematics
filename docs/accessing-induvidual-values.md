## Accessing Individual Values

To evaluate kinematic quantities at a specific value, use `at_value()`.

This method automatically handles multi-valued solutions and always returns lists.

### Syntax

```python
from reaction_kinematics import Reaction

rxn = Reaction("p", "3H", "n", "3He")
# Proton + Tritium Reaction

r = rxn.at_value("theta3", 2.13, ek=1.2)
# 2.13 is in radians
print(r)

```
This will return:
```python
{'coscm': [-0.9813047694699246], 
'theta_cm': [2.9479854759117416], 
'theta3': [2.13], 
'theta4': [0.054070851022814864], 
'e3': [0.005237930117122223], 
'e4': [0.43100601988343396], 
'v3': [0.003338681989975168], 
'v4': [0.017517730726598162], 
'p3': [3.136927647029636], 
'p4': [49.204196838704895]}
```
Parameters:

* `x_name` : Independent variable (e.g. `"theta3"`, `"theta4"`, `"theta_cm"`, `"coscm"`)
* `x`      : Value at which to evaluate (radians for angles)
* `ek`     : Beam kinetic energy in MeV (required)
* `y_names`: Dependent variables to return — string, list, or `None` for all
