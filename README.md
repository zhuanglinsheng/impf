# impf

This is a Python package for solving implicit functions, and more.



## Features

- This is essentially a pure C90 library bound to Python.
- No third-party dependencies.
- C APIs are designed in the explicit way, without using `typedef`.
- The Python component handles the UI work.

Currently, you can use this package for

- Linear Programming
- Interpolation



## Build & Install

To build and install this package from source, enter the root folder and type

```shell
pip3 install -e . -v
```

Then, the installation should be done.



## Example

### 1. Linear Programming

```python
from impf import optm

obj = optm.LinearObjective([3, 2], 'max')
constraint_1 = optm.LinearConstraint([1, 1], '<=', 9)
constraint_2 = optm.LinearConstraint([3, 1], '<=', 18)
constraint_3 = optm.LinearConstraint([1, 0], '<=', 7)
constraint_4 = optm.LinearConstraint([0, 1], '<=', 6)
constraints = [constraint_1, constraint_2, constraint_3, constraint_4]

prob = optm.LinearProgramming(obj, constraints)
prob
```

<div>
max&nbsp;&nbsp;[3.0, 2.0]<br>
s.t.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[1.0, 1.0] <= 9.0<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.0, 1.0] <= 18.0<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[1.0, 0.0] <= 7.0<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[0.0, 1.0] <= 6.0<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;all variables >= 0<br>
</div>

```python
prob.solve()
```

<div>
{'x': [4.5, 4.499999999999999], 'state': 'Success'}
</div>

