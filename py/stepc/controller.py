# STEPC
# A fast and extensible linear system and controller simulator
# for evaluation of predictive controllers
#
# Jon Sowman 2014 <j.sowman@soton.ac.uk>

import numpy as np
import cvxopt

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
        self.__Hu = 0
        self.__target = 0

        self.__G = 0

        # The problem Hessian matrix
        self.__H = 0

        # Store the previous input
        self.u_last = 0

    def set_prediction_horizon(self, Hp):
        self.__Hp = Hp

    def set_control_horizon(self, Hu):
        self.__Hu = Hu

    def set_target(self, target):
        self.__target = target

    def reset(self):
        self.u_last = 0

    def generate_matrices(self):
        """
        Use a linear MPC technique along with the weights P, Q and R
        to construct and solve a linear MPC problem to give the
        required control move.

        The problem formulation in this method is based on Maciejowski (2002)
        and is detailed in pages 74-76.

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

        # Rbar = diag(R, R, ..., R) with dimension lHu x lHu (l = numinputs)
        Rbar = np.zeros([self.__Hu * self.__sys.numinputs, self.__Hu *
                        self.__sys.numinputs])
        for idx in range(self.__Hu):
            start = idx * self.__sys.numinputs
            end = start + self.__sys.numinputs
            Rbar[start:end, start:end] = R
        
        # Psi = [A; A^2; A^3 ... ] with dimension mHp x n
        psi = np.zeros([self.__sys.order * self.__Hp, self.__sys.order])
        for idx in range(self.__Hp):
            start = idx * self.__sys.order
            end = start + self.__sys.order
            # Each block is A ^ (index+1) (since index starts at 0)
            psi[start:end, 0:self.__sys.order] \
                = np.linalg.matrix_power(self.__sys.A, idx+1)

        # gamma = [B; ...; sum_{i=0}^{Hp-1} A^i * B]
        gamma = self.__sys.B.copy()
        accum = self.__sys.B.copy()
        for idx in range(1, self.__Hp):
            exp_A = np.linalg.matrix_power(self.__sys.A, idx)
            accum += exp_A.dot(self.__sys.B)
            gamma = np.vstack((gamma, accum))

        # phi = gamma in first col, then shifted down by one for each column
        # until a total of Hu columns exists
        phi = np.zeros([self.__sys.order * self.__Hp,
                        self.__sys.numinputs * self.__Hu])
        for idx in range(self.__Hu):
            rowstart = idx * self.__sys.order
            rowend = (self.__Hp - idx) * self.__sys.order
            colstart = idx * self.__sys.numinputs
            colend = colstart + self.__sys.numinputs
            phi[rowstart:self.__sys.order * self.__Hp, colstart:colend] \
                    = gamma[0:rowend, :]

        # tau = [target; target; ... ] for all of Hp
        tau = self.__target.copy()
        for idx in range(1, self.__Hp):
            tau = np.vstack((tau, self.__target))

        # Store psi, gamma and phi
        self.__psi = psi
        self.__phi = phi
        self.__gamma = gamma
        self.__tau = tau

        # G = 3 theta.T Qbar
        self.__G = 2 * phi.T.dot(Qbar)

        # H = 2(Rbar + B.T Qbar B)
        self.__H = phi.T.dot(Qbar).dot(phi)
        self.__H = Rbar + self.__H
        self.__H = self.__H * 2

    def controlmove(self, x0):
        """
        Use Linear MPC to calculate a control move, returning the first input.

        This method uses the cvxopt quadratic programme (QP) solver (which uses
        a primal-dual interior point method) to find the optimal system input
        given the internal system model (self.__sys). The optimisation finds an
        input vector over the control horizon (self.__Hu) but only the first
        input is returned and applied to the plant.
        """

        # Find the 'tracking error' of the controller
        # error = tau - psi*x0 - gamma*u_old
        error = self.__tau - self.__psi.dot(x0)
        error -= self.__gamma.dot(self.u_last)

        # Let Gx be G.dot(error)
        Geps = self.__G.dot(error)
        
        # Need to convert to 'cvxopt' matrices instead of np arrays
        cvx_H = cvxopt.matrix(self.__H)
        cvx_Geps = cvxopt.matrix(Geps)

        # Run the optimiser (note the negative here, see Maciejowski eq. 3.10)
        results = cvxopt.solvers.qp(cvx_H, -cvx_Geps)

        # Extract result and turn it back to an np array
        uvect = np.array(results['x'])

        # Return u0, the first control input and store it for next iteration
        u0 = np.array([uvect[0]]) + self.u_last
        self.u_last = u0
        return u0
