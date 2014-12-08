# STEPC
# A fast and extensible linear system and controller simulator
# for evaluation of predictive controllers
#
# Jon Sowman 2014 <j.sowman@soton.ac.uk>

from ..stepc.solver import *

import numpy as np
from pprint import pprint


def test_forwards_euler_single():
    """
    Create an initial condition and then use forwards Euler
    to integrate it forwards
    """
    solver = ForwardsEulerSolver()
    x0 = np.array([[0], [1], [2]])
    xdot = np.array([[-1], [1], [0]])
    x = solver.step(x0, xdot, 1)

    assert np.allclose(x, np.array([[-1], [2], [2]]))


def test_forwads_euler_multiple():
    """
    Create an initial condition and then use forwards Euler
    to integrate it forwards, using 1000 time steps to test numerical
    stability
    """
    solver = ForwardsEulerSolver()
    x = np.array([[0], [1], [2]])
    xdot = np.array([[-0.1], [0.01], [0.05]])
    for i in range(1000):
        x = solver.step(x, xdot, 1)

    assert np.allclose(x, np.array([[-100.], [11.], [52.]]))
