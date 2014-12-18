# STEPC
# A fast and extensible linear system and controller simulator
# for evaluation of predictive controllers
#
# Jon Sowman 2014 <j.sowman@soton.ac.uk>

import numpy as np


class LinSystem(object):
    def __init__(self, order, numinputs, numoutputs):
        self.A = np.empty([order, order])
        self.B = np.empty([order, numinputs])
        self.C = np.empty([numoutputs, order])
        self.order = order
        self.numinputs = numinputs
        self.numoutputs = numoutputs

    order = 0
    numinputs = 0
    numoutputs = 0
