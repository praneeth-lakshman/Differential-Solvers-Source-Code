import sys, os
script_dir = os.path.dirname(os.path.abspath('__file__'))
module_dir = os.path.join(script_dir, "build")
sys.path.append(module_dir)
import ODESolvers as ode
import math

