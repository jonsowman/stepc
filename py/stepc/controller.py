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

    def set_plant(self, sys):
        self.__sys = sys

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

    def generate_constraints(self, ulim, dulim, zlim):
        """
        Generate the linear inequality constraint matrices Ax < b
        for the quadratic programming problem. Constraints are to be phrased as
        constraints on delta-u (i.e. changes in the control inputs).

        This formulation is detailed in Maciejowski pp. 81-83.

        ulim contains 2 colums and as many rows as there are plant inputs. The
        two columns are the lower and upper limits on u respectively. The same
        applies for dulim (delta-u limits) and zlim (controlled variable
        limits).

        """

        # Make sure generate_matrices has been run as we rely on some of its
        # results
        try:
            self.__phi
        except NameError:
            print "You need to run generate_matrices() before constraints \
                    can be added"

        # Check that limits have 2 columns each
        assert (ulim.shape[1] == 4), "ulim must have exactly 4 columns, \
                has %d" % ulim.shape[1]
        assert (dulim.shape[1] == 4), "dulim must have exactly 4 columns, \
                has %d" % dulim.shape[1]
        assert (zlim.shape[1] == 4), "zlim must have exactly 4 columns, \
                has %d" % zlim.shape[1]

        # Extract the upper and lower limits, their weights, and force them to
        # be 2D row vectors
        ulow = ulim[:, 0].reshape(1, -1)
        ulow_weight = ulim[:, 1].reshape(1, -1)
        uhigh = ulim[:, 2].reshape(1, -1)
        uhigh_weight = ulim[:, 3].reshape(1, -1)

        dulow = dulim[:, 0].reshape(1, -1)
        dulow_weight = dulim[:, 1].reshape(1, -1)
        duhigh = dulim[:, 2].reshape(1, -1)
        duhigh_weight = dulim[:, 3].reshape(1, -1)

        zlow = zlim[:, 0].reshape(1, -1)
        zlow_weight = zlim[:, 1].reshape(1, -1)
        zhigh = zlim[:, 2].reshape(1, -1)
        zhigh_weight = zlim[:, 3].reshape(1, -1)

        # Verify that the constraint vectors given are of the right dimension
        assert (np.size(dulow) == self.__sys.numinputs), "Size of dulow \
                 constraint vector (%d) and number of system inputs \
                 (%d) must be the same" % (np.size(dulow),
                                           self.__sys.numinputs)
        assert (np.size(duhigh) == self.__sys.numinputs), "Size of duhigh \
                constraint vector (%d) and number of system inputs \
                (%d) must be the same" % (np.size(duhigh),
                                          self.__sys.numinputs)
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

        # SOFT CONSTRAINTS
        # Stack weights in the same order as the constraint coefficients
        weights = np.vstack((np.tile(uhigh_weight.T, (self.__Hu, 1)),
                             np.tile(ulow_weight.T, (self.__Hu, 1)),
                             np.tile(duhigh_weight.T, (self.__Hu, 1)),
                             np.tile(dulow_weight.T, (self.__Hu, 1)),
                             np.tile(zhigh_weight.T, (self.__Hp, 1)),
                             np.tile(zlow_weight.T, (self.__Hp, 1))))

        # Check for zero-weights since they break the QP and should be removed
        assert weights.all(), "All weights must be positive real numbers"

        # Store
        self.weights = weights

    def soften_constraints(self, H, G, A, B, weights):
        """
        The QP objective function is min_{x} x'Hx + Gx.

        We are given a set of inequality constraints Ax < b, and a vector of
        weights that are associated with these constraints.

        We should adjust the phrasing of the inequality constraints such that
        they are 'softened' as specified by the 'weights' vector. Weights of 0
        are to be discarded (the constraint is "infinitely soft" and need not
        be respected at all. Weights of +inf are "hard" and are never violated
        by the QP. Weights of < 0 are not permitted and should fatally error.
        """

        # Some sanity checking
        assert weights.all(), "All weights must be positive real numbers"
        assert not np.isnan(weights).any(), "All weights must be numerical"
        assert not np.isneginf(weights).any(), "Negative infinity is an invalid \
                weight"

        # Find the indices of the soft and hard constraints
        soft_idx = np.where(np.logical_and(weights != 0, weights != np.inf))[0]
        num_soft = np.size(soft_idx)

        # Get the soft constraint weights
        weights_soft_only = np.take(weights, soft_idx)

        # Number of constraints
        (num_constraints, num_variables) = A.shape
        minusones = np.zeros([num_constraints, num_soft])
        for i in range(num_soft):
            minusones[soft_idx[i], i] = -1

        # Only bother doing this if there are any soft constraints
        if num_soft > 0:
            # Add rows to A and B for the soft constraints
            A = np.vstack((np.hstack((A, minusones)),
                           np.hstack((np.zeros([num_soft, num_variables]),
                                     -np.eye(num_soft)))))
            B = np.vstack((B, np.zeros([num_soft, 1])))

        if weights_soft_only.size:
            # Determine the size of the current objective function Hessian
            [Hrows, Hcols] = H.shape
            # Append the slack variables (epsilon) to the objective function
            H = np.vstack((np.hstack((H, np.zeros([Hrows, num_soft]))),
                           np.hstack((np.zeros([num_soft, Hcols]),
                                      np.diagflat(weights_soft_only)))))

            # Add zeros to G (this is a 2-norm subtlety - change for 1-norm)
            G = np.vstack((G, np.zeros([num_soft, 1])))

        return (H, G, A, B)

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

        # Soften the constraints
        (H_soft, G_soft, A_soft, B_soft) = \
            self.soften_constraints(self.__H,
                                    Geps, constraints_lhs_stacked,
                                    constraints_rhs_stacked,
                                    self.weights)

        # Need to convert to 'cvxopt' matrices instead of np arrays
        cvx_H = cvxopt.matrix(H_soft)
        cvx_Geps = cvxopt.matrix(G_soft)
        cvx_constraints_lhs_soft = cvxopt.matrix(A_soft)
        cvx_constraints_rhs_soft = cvxopt.matrix(B_soft)

        # Run the optimiser (note the negative here, see Maciejowski eq. 3.10)
        cvxopt.solvers.options['show_progress'] = False
        results = cvxopt.solvers.qp(cvx_H, -cvx_Geps,
                                    cvx_constraints_lhs_soft,
                                    cvx_constraints_rhs_soft)

        # Extract result and turn it back to an np array
        uvect = np.array(results['x'])

        # Return u0, the first control input and store it for next iteration
        u0 = np.array([uvect[0]]) + self.u_last
        self.u_last = u0
        return u0
