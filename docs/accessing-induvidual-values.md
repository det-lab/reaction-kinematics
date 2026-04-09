## Accessing Individual Values

To evaluate kinematic quantities at a specific value, use `at_value()`.

This method automatically handles multi-valued solutions and always returns lists.

### Syntax

```python
rxn.at_value(x_name, x_value, y_names=None)
```

Parameters:

* `x_name` : Independent variable (e.g. `"theta4"`, `"theta_cm"`, `"coscm"`)
* `x_value`: Value at which to evaluate
* `y_names`: Dependent variables (string or list)