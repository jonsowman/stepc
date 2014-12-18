# STEPC
# A fast and extensible linear system and controller simulator
# for evaluation of predictive controllers
#
# Jon Sowman 2014 <j.sowman@soton.ac.uk>

import numpy as np


class Solver(object):
    def step(self):
        return 0


class ForwardsEulerSolver(Solver):
    def step(self, x0, xdot, Ts):
        return x0 + np.dot(Ts, xdot)
