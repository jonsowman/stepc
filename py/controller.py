# STEPC
# A fast and extensible linear system and controller simulator
# for evaluation of predictive controllers
#
# Jon Sowman 2014 <j.sowman@soton.ac.uk>

import numpy as np

from linsystem import LinSystem


class Controller:
    def __init__(self):
        self.sys = 0

    def controlmove(self):
        """
        Returns a null control move. This should be overridden by
        derived classes such as the PIDController
        """
        return np.empty([self.sys.numinputs, 1])


class PIDController(Controller):
    def __init__(self):
        self.integrator = 0
        self.old_error = 0
        self.kp = 0
        self.ki = 0
        self.kd = 0

    def controlmove(self, y):
        """
        Calculate a control move based on the system output vector y
        """

        # Find the current error
        e = self.target - y

        # Accumulate the error over time
        self.integrator += e

        # Find the control action in PID form
        u = np.dot(e, self.kp)
        u += np.dot(self.integrator, self.ki)
        u += np.dot((e - self.old_error), self.kd)

        # Store the old error
        self.old_error = e

        return u

    def set_kp(self, kp):
        self.kp = kp

    def set_ki(self, ki):
        self.ki = ki

    def set_kd(self, kd):
        self.kd = kd

    def set_target(self, target):
        self.target = target
