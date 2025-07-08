# Differential Solvers Source Code

This is a package to numerically solve 1st Order ODEs with initial values entirely written by me in C++, with a Python wrapper so it can be easily used in your code repos. Everything is implemented in C++17, except for the Forward/Backward/Central finite difference methods, which are in native Python. This provides code that is marginally faster and more accurate than the default SciPy solver, but the difference is neglible at best.

Check out the cpp_files subdirectory to find the original C++ files from which I added to the main.cpp file (if you just want to study one of the ODE methods)

## What's Included

In order of least accurate to most accurate

- Euler's Method
- Backwards Euler Method
- Runge-Kutta 4th Order
- Runge-Kutta-Fehlberg (RK4(5))
- Dormand Prince

I intend to add the theory behind the equations, either as a video series or a Jupyter Notebook.

## How to Install (Only for Linux systems)

Ensure that your Python version is 3.10 - if you are on a different version, please [recompile the binary!](#Recompile) Also, ensure that numpy is installed (pip install numpy)

1. Open your terminal in a folder containing your existing code and run this command:

```shell
git clone https://github.com/praneeth-lakshman/differential_solvers.git .
```

2. In your Python code, import the file as so:

```python
import ODE as ode
```

That's it! Now you can use the module with all of its methods and classes

## Manual

(ODE is ode, as it assuming you are following the install instructions)
| ODE Method | Implementation |
| ------- | -------------- |
| Euler's Method | ode.Explicit.Euler(y, t_span, h, dydt) |
| Backwards Euler Method | ode.Implicit.BackwardsEuler(y, t_span, h, dydt) |
| Runge-Kutta 4th Order | ode.Explicit.RK4(y, t_span, h, dydt) |
| Runge-Kutta-Fehlberg | ode.Explicit.RKF(y, t_span, h, dydt) |
| Dormand-Prince | ode.Explicit.DormandPrince(y, t_span, h, dydt) |

All of the following (except forward finite difference) is implemented in Python as it is very computationally inexpensive
| Miscellaneous Method | Implementation |
| ------------------------------- | -------------------------------------------------- |
| Forward finite difference | ode.ODEMethods.finite_difference(f, t, dt) |
| Backward finite difference | ode.ODEMethods.backward_difference(f, t, dt) |
| Central finite difference | ode.ODEMethods.central_difference(f, t, dt) |
| Multivariable finite difference | ode.ODEMethods.finite_difference_of_y(f, y, t, dt) |

For the last one, "y" is the variable you want to differentiate w.r.t and "t" is a variable input to f that we are NOT differentiating.

# Recompile

Ensure CMake is installed, and other necessary tools

1. Remove the build directory, and move into it

```shell
rm -rf build/* && cd build;
```

2. CMake home directory

```shell
cmake ..
```

3. Make the .so file

```shell
make
```

4. Note the name of the .so file generated

5. Go to ODE.py, and change Line 4 to that name
