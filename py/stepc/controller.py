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

    def generate_constraints(self, ulow, uhigh, dulow, duhigh, zlow, zhigh):
        """
        Generate the linear inequality constraint matrices Ax < b
        for the quadratic programming problem. Constraints are to be phrased as
        constraints on delta-u (i.e. changes in the control inputs).

        This formulation is detailed in Maciejowski pp. 81-83.

        ulow is a vector specifying the lower box constraint for each control
        input u (i.e. this vector has length equal to the number of controlled
        inputs of the plant). A similar situation applies for uhigh.

        dulow and duhigh are the same thing but for changes in the input that
        are permitted between one time step and the next.
        """

        # Make sure generate_matrices has been run as we rely on some of its
        # results
        try:
            self.__phi
        except NameError:
            print "You need to run generate_matrices() before constraints \
                    can be added"

        # Verify that the constraint vectors given are of the right dimension
        assert (np.size(dulow) == self.__sys.numinputs), "Size of dulow \
                 constraint vector (%d) and number of system inputs \
                 (%d) must be the same" % (np.size(dulow), self.__sys.numinputs)
        assert (np.size(duhigh) == self.__sys.numinputs), "Size of duhigh \
                constraint vector (%d) and number of system inputs \
                (%d) must be the same" % (np.size(duhigh), self.__sys.numinputs)
        assert (np.size(ulow) == self.__sys.numinputs), "Size of ulow \
               constraint vector (%d) and number of system inputs \
               (%d) must be the same" % (np.size(ulow), self.__sys.numinputs)
        assert (np.size(uhigh) == self.__sys.numinputs), "Size of uhigh \
                constraint vector (%d) and number of system inputs \
                (%d) must be the same" % (np.size(uhigh), self.__sys.numinputs)
        assert (np.size(zlow) == self.__sys.order), "Size of zlow \
               constraint vector (%d) and order of plant \
               (%d) must be the same" % (np.size(zlow), self.__sys.order)
        assert (np.size(zhigh) == self.__sys.order), "Size of zhigh \
                constraint vector (%d) and oder of plant \
                (%d) must be the same" % (np.size(zhigh), self.__sys.order)

        # These matrices are used during the generation of more than one set of
        # constraint matrices
        
        # I matrices in A are m x m (m = number of controlled inputs)
        I = np.eye(self.__sys.numinputs)
        block = np.vstack((I, -I))

        # Create a zero matrix of the same size as 'block'
        z = np.zeros([self.__sys.numinputs * 2, self.__sys.numinputs])

        # DELTA-INPUT CONSTRAINTS
        # Start with an empty E
        self.W = np.zeros([self.__sys.numinputs * self.__Hu * 2, 0])
        # Create the full E matrix recursively
        for col_idx in range(self.__Hu):
            col_top = np.tile(z, (col_idx, 1))
            col_bot = np.tile(z, (self.__Hu - col_idx - 1, 1))
            col = np.vstack((col_top, block, col_bot))
            self.W = np.hstack((self.W, col))

        # Generate w (little-w) which is the RHS of the delta-u inequality
        w_single = np.vstack((duhigh.T, -dulow.T))
        self.w = np.tile(w_single, (self.__Hu, 1))

        # INPUT CONSTRAINTS
        # Start with an empty F
        self.F = np.zeros([self.__sys.numinputs * self.__Hu * 2, 0])
        # Create the full F matrix recursively
        for col_idx in range(self.__Hu):
            col_top = np.tile(z, (col_idx, 1))
            col_bot = np.tile(block, (self.__Hu - col_idx, 1))
            col = np.vstack((col_top, col_bot))
            self.F = np.hstack((self.F, col))

        # Extract F1 since we need it for the RHS (3.37 in Maciejowski)
        self.F1 = self.F[:, 0:self.__sys.numinputs]

        # Now generate 'f' (little-f)
        f_single = np.vstack((-uhigh.T, ulow.T))
        self.f = np.tile(f_single, (self.__Hu, 1))

        # STATE CONSTRAINTS
        # Gamma is similar to W in the above, but could be of different
        # dimension since the number of controller outputs (Z) is not
        # necessarily the same as the number of plant inputs (U)
        # First, regenerate the construction matrices
        I = np.eye(self.__sys.order)
        block = np.vstack((I, -I))
        z = np.zeros([self.__sys.order * 2, self.__sys.order])

        self.Gamma = np.zeros([self.__sys.order * self.__Hp * 2, 0])
        # Create the full Gamma matrix recursively
        for col_idx in range(self.__Hp):
            col_top = np.tile(z, (col_idx, 1))
            col_bot = np.tile(z, (self.__Hp - col_idx - 1, 1))
            col = np.vstack((col_top, block, col_bot))
            self.Gamma = np.hstack((self.Gamma, col))

        # Now construct g (little-g) which forms part of the RHS of the state
        # (controlled output) inequality
        g_single = np.vstack((-zhigh.T, zlow.T))
        self.g = np.tile(g_single, (self.__Hp, 1))

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

        # Sort the constraints
        # Input constraints
        u_rhs = -self.F1.dot(self.u_last) - self.f

        # State constraints
        z_lhs = self.Gamma.dot(self.__phi)
        z_rhs = -self.Gamma.dot(self.__psi.dot(x0) -
                self.__gamma.dot(self.u_last)) - self.g

        # Stack the constraints
        constraints_lhs_stacked = np.vstack((self.F, self.W, z_lhs))
        constraints_rhs_stacked = np.vstack((u_rhs, self.w, z_rhs))

        # Need to convert to 'cvxopt' matrices instead of np arrays
        cvx_H = cvxopt.matrix(self.__H)
        cvx_Geps = cvxopt.matrix(Geps)
        cvx_constraints_lhs_stacked = cvxopt.matrix(constraints_lhs_stacked)
        cvx_constraints_rhs_stacked = cvxopt.matrix(constraints_rhs_stacked)

        # Run the optimiser (note the negative here, see Maciejowski eq. 3.10)
        cvxopt.solvers.options['show_progress'] = False
        results = cvxopt.solvers.qp(cvx_H, -cvx_Geps,
                cvx_constraints_lhs_stacked, cvx_constraints_rhs_stacked)

        # Extract result and turn it back to an np array
        uvect = np.array(results['x'])

        # Return u0, the first control input and store it for next iteration
        u0 = np.array([uvect[0]]) + self.u_last
        self.u_last = u0
        return u0
