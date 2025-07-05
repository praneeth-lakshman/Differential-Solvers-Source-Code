import sys, os
from collections.abc import Callable
script_dir = os.path.dirname(os.path.abspath('__file__'))
module_dir = os.path.join(script_dir, "build")
sys.path.append(module_dir)
import ODESolvers as ode
import math
import numpy as np

class Explicit:
    def Euler(y: float, t_span: (float, float), h: float, dydt):
        """Use the (forward) Euler method to approximate a solution to IVP. Returns list of y values

        Args:
            y (float): initial value of y
            t_span (float, float): Tuple of range (inclusive at both ends) of t
            h (float): step size
            dydt (Callable): The equation of the ODE
        """
        t_values = np.linspace(t_span[0], t_span[1], num=math.floor((t_span[1]-t_span[0])/h)+1)
        y_vals = []
        for t in t_values:
            y_vals.append(y)
            y = ode.Euler().next_step(y, h, dydt, t)
        return y_vals

def func(y, t):
    return t + 3*y

print(Explicit.Euler(1, (0,3), 0.1, func))
        
