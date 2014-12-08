# STEPC
# A fast and extensible linear system and controller simulator
# for evaluation of predictive controllers
#
# Jon Sowman 2014 <j.sowman@soton.ac.uk>

from ..stepc.controller import LinearMPCController
from ..stepc.linsystem import LinSystem

import numpy as np


class TestLinearMPC:

    def setup(self):
        """
        Create all the bits needed to test the ability of a LinearMPCController
        to derive and formulate a solution framework, and then find a solution
        to a toy problem.
        """

        self.__sys = LinSystem(2, 1, 1)
        self.__sys.A = np.eye(2)
        self.__sys.B = np.array([[1], [2]])

    def test_generate_F(self):
        mpc = LinearMPCController(self.__sys)
        mpc.P = mpc.Q = mpc.R = 1
        mpc.set_prediction_horizon(1)
        mpc.generate_matrices()

        expected = np.array([[2, 4]])

        assert np.array_equal(expected, mpc._LinearMPCController__F)
