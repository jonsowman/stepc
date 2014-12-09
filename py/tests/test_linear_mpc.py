# STEPC
# A fast and extensible linear system and controller simulator
# for evaluation of predictive controllers
#
# Jon Sowman 2014 <j.sowman@soton.ac.uk>

from ..stepc.controller import LinearMPCController
from ..stepc.linsystem import LinSystem

import numpy as np


class TestLinearMPC(object):

    def setup(self):
        """
        Create all the bits needed to test the ability of a LinearMPCController
        to derive and formulate a solution framework, and then find a solution
        to a toy problem.
        """

        self.__sys = LinSystem(2, 1, 1)
        # Test system is from 12mr (mass/spring)
        self.__sys.A = np.array([[0.1, 1], [1, -0.1]])
        self.__sys.B = np.array([[0], [0.1]])
    setup.__test__ = False

    def teardown(self):
        pass
    teardown.__test__ = False

    def test_generate_F(self):
        mpc = LinearMPCController(self.__sys)
        mpc.P = mpc.Q = mpc.R = 1
        mpc.set_prediction_horizon(1)
        mpc.generate_matrices()

        expected_F = np.array([[0.2, -0.02]])

        assert np.allclose(expected_F, mpc._LinearMPCController__F)

    def test_generate_G(self):
        mpc = LinearMPCController(self.__sys)
        mpc.P = mpc.Q = mpc.R = 1
        mpc.set_prediction_horizon(1)
        mpc.generate_matrices()

        expected_G = np.array([[2.02]])

        assert np.allclose(expected_G, mpc._LinearMPCController__G)
