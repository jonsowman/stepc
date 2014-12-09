# STEPC
# A fast and extensible linear system and controller simulator
# for evaluation of predictive controllers
#
# Jon Sowman 2014 <j.sowman@soton.ac.uk>

import numpy as np

from linsystem import LinSystem
from pprint import pprint


class Controller(object):
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


class LinearMPCController(Controller):
    def __init__(self, sys):
        self.P = 0
        self.Q = 0
        self.R = 0
        self.__sys = sys
        self.__Hp = 0

        self.__F = 0
        self.__G = 0

    def set_prediction_horizon(self, Hp):
        self.__Hp = Hp

    def generate_matrices(self):
        """
        Use a linear MPC technique along with the weights P, Q and R
        to construct and solve a linear MPC problem to give the
        required control move.
        """

        # Assert some things before we break maths
        if not self.P or not self.Q or not self.R:
            raise AssertionError("MPC Matrices have not been set")

        if not self.__sys:
            raise AssertionError("A model has not been attached to \
                                the MPC controller")

        if not self.__Hp:
            raise AssertionError("A prediction horizon has not been set")

        # Q & P diagonal with number of elements = number of states
        Q = self.Q * np.eye(self.__sys.order)
        P = self.P * np.eye(self.__sys.order)
        # Same for R
        R = self.R * np.eye(self.__sys.numinputs)

        # Qbar = diag(Q, Q, Q, ..., P)
        Qbar = np.zeros([self.__Hp * self.__sys.order, self.__Hp *
                        self.__sys.order])
        for idx in range(self.__Hp):
            start = idx * self.__sys.order
            end = start + self.__sys.order
            Qbar[start:end, start:end] = Q

        # Add P to the final diagnonal element of Qbar
        Qbar[-self.__Hp:-1, -self.__Hp:-1] = self.P

        # Rbar = diag(R, R, ..., R)
        Rbar = np.zeros([self.__Hp * self.__sys.numinputs, self.__Hp *
                        self.__sys.numinputs])
        for idx in range(self.__Hp):
            start = idx * self.__Hp
            end = start + self.__sys.numinputs
            Rbar[start:end, start:end] = R
        
        # Abar = [A A^2 A^3 ... ]
        Abar = np.zeros([self.__sys.order * self.__Hp, self.__sys.order])
        for idx in range(self.__Hp):
            start = idx * self.__sys.order
            end = start + self.__sys.order
            # Each block is A ^ (index+1) (since index starts at 0)
            Abar[start:end, 0:self.__sys.order] \
                = np.linalg.matrix_power(self.__sys.A, idx+1)

        # Bbar = [A, 0; AB, A] for example (Hp=2)
        Bbar = np.zeros([self.__sys.order * self.__Hp,
                        self.__sys.numinputs * self.__Hp])

        # log_ = "logical", i.e. the submatrix blocks
        # phy_ = "physical", i.e. the way the matrix is laid
        #         out in memory
        phy_row = phy_col = log_row = log_col = 0
        for log_row in range(self.__Hp):
            for log_col in range(self.__Hp):
                phy_row += self.__sys.order
                phy_col += self.__sys.numinputs

                if(log_row >= log_col):
                    phy_row_start = log_row * self.__sys.order
                    phy_row_end = phy_row_start + self.__sys.order
                    phy_col_start = log_col * self.__sys.numinputs
                    phy_col_end = phy_col_start + self.__sys.numinputs
                    Bbar[phy_row_start:phy_row_end, phy_col_start:phy_col_end] \
                        = np.linalg.matrix_power(self.__sys.A,
                                log_row - log_col).dot(self.__sys.B)

        # F = 2 B.T Qbar A
        self.__F = 2 * Bbar.T.dot(Qbar).dot(Abar)

        # G = 2(Rbar + B.T Qbar B)
        self.__G = Bbar.T.dot(Qbar).dot(Bbar)
        self.__G = Rbar + self.__G
        self.__G = self.__G * 2
