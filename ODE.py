import sys, os
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
    def RK4(y: float, t_span: (float, float), h: float, dydt):
        """Use the Runge Kutta 4th order method to approximate a solution to IVP. Returns list of y values

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
            y = ode.RK4().next_step(y, h, dydt, t)
        return y_vals
    def RKF(y: float, t_span: (float, float), h: float, dydt):
        """Use the Runge-Kutta-Fehlberg method to approximate a solution to IVP. Returns tuple of lists: (t_vals, y_vals)

        Args:
            y (float): initial value of y
            t_span (float, float): Tuple of range (inclusive at both ends) of t
            h (float): step size
            dydt (Callable): The equation of the ODE
        """
        t = t_span[0]
        t_vals = []
        y_vals = []
        while t <= t_span[1]:
            t_vals.append(t)
            y_vals.append(y)
            result = ode.RKF().next_step(y, h, dydt, t)
            if result.h < 1e-6:
                break
            if result.accepted:
                y = result.y
                t += h
            h = result.h
        return (t_vals, y_vals)

    def DormandPrince(y: float, t_span: (float, float), h: float, dydt):
        """Use the Dormand-Prince method to approximate a solution to IVP. Returns tuple of lists: (t_vals, y_vals)

        Args:
            y (float): initial value of y
            t_span (float, float): Tuple of range (inclusive at both ends) of t
            h (float): step size
            dydt (Callable): The equation of the ODE
        """
        t = t_span[0]
        t_vals = []
        y_vals = []
        while t <= t_span[1]:
            t_vals.append(t)
            y_vals.append(y)
            result = ode.DormandPrince().next_step(y, h, dydt, t)
            if result.h < 1e-6:
                break
            if result.accepted:
                y = result.y
                t += h
            h = result.h
        return (t_vals, y_vals)

class Implicit:
    def BackwardsEuler(y: float, t_span: (float, float), h: float, dydt):
        """Use the (backwards) Euler method to approximate a solution to IVP. Returns list of y values

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
            y = ode.BackwardsEuler().next_step(y, h, dydt, t)
        return y_vals

class ODEMethods:
    def finite_difference_of_y(f, y, t, dt=0.001):
        """Forward finite difference of a multivariable function f that depends on y and t, to calculate df/dy

        Args:
            f (function): a function that depends on y and t, from which df/dy is needed
            y (float): y value at point (differentiation is wrt this variable)
            t (float): t value at point (differentiation is NOT wrt this variable)
            dt (float, optional): time step. Defaults to 0.001.

        Returns:
            float: value of function f at next time step
        """
        return ode.finite_difference(y, f, t, dt)
    def finite_difference(f, t, dt=0.001):
        """Forward finite difference of a function f that depends on t ONLY, to calculate df/dt

        Args:
            f (function): a function that depends on t only, from which df/dt is required
            t (float): t value at point
            dt (float, optional): time step. Defaults to 0.001.

        Returns:
            float: value of function f at the next time step 
        """
        return (f(t+dt) - f(t)) / dt
    def find_zero(self, f, t, y, e=1e-6):
        """Implements Newton-Raphson method for a function f(y, t)

        Args:
            f (function): a function that depends on y and t, that equals 0
            t (float): t value at point
            y (float): y value at point
            e (float, optional): Sensitivity of root finding algorithm. The lower this value, the more iterations. Defaults to 1e-6.
        """
        return ode.find_zero(f, t, y, e)
    def backward_difference(f, y, dt=0.001):
        """Backward difference of a function f that depends on t ONLY, to calculate df/dt

        Args:
            f (function): a function that depends on t only, from which df/dt is required
            t (float): t value at point
            dt (float, optional): time step. Defaults to 0.001.

        Returns:
            float: value of function f at the next time step 
        """
        return (f(t) - f(t-dt))/dt
    def central_difference(f, y, dt=0.002):
        """Central finite difference of a function f that depends on t ONLY, to calculate df/dt

        Args:
            f (function): a function that depends on t only, from which df/dt is required
            t (float): t value at point
            dt (float, optional): Total time step (this means that this method goes dt/2 forwards and dt/2 backwards). Defaults to 0.001.

        Returns:
            float: value of function f at the next time step 
        """
        return (f(t+dt/2) - f(t-dt/2))/dt
    