# STEPC
# A fast and extensible linear system and controller simulator
# for evaluation of predictive controllers
#
# Jon Sowman 2014 <j.sowman@soton.ac.uk>

import numpy as np
import scipy.linalg
from linsystem import LinSystem


def c2d(sys, Ts):
    """
    Convert a continuous time system described by
        xdot = Ax + Bu
    into its discrete time equivalent. The sample time used for this conversion
    is given by Ts.

    Based on lecture notes from Dr D. S. Laila, University of Southampton
    """

    # Extract data from the LinSystem
    A = sys.A
    B = sys.B

    [rows_a, cols_a] = A.shape
    [rows_b, cols_b] = B.shape

    augmented = np.vstack((np.hstack((A, B))*Ts,
                           np.zeros([cols_b, cols_a + cols_b])))
    result = scipy.linalg.expm(augmented)

    phi = result[0:cols_a, 0:cols_a]
    gamma = result[0:cols_a, cols_a:cols_a+cols_b]

    # Construct and return a new LinSystem
    sysd = LinSystem(sys.order, sys.numinputs, sys.numoutputs)
    sysd.A = phi
    sysd.B = gamma
    sysd.C = sys.C

    return sysd
