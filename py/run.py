# STEPC
# A fast and extensible linear system and controller simulator
# for evaluation of predictive controllers
#
# Jon Sowman 2014 <j.sowman@soton.ac.uk>

import numpy as np

from linsystem import LinSystem
from simulator import Simulator

# Create the system
sys = LinSystem(2, 1, 1)
sys.A[0, 0] = 0
sys.A[0, 1] = 1
sys.A[1, 0] = -1
sys.A[1, 1] = 0

sys.B[0, 0] = 0
sys.B[1, 0] = -1

sys.C[0, 0] = 1
sys.C[0, 1] = 0

# Make a simulator and set parms
sim = Simulator()
sim.timestep = 0.01
sim.endtime = 10

# Initial condition
x0 = np.matrix([[-1], [0]])

# Go for it
x = sim.simulate(sys, x0)

# Results
print "Final state is [%.5f, %.5f]" % (x[0], x[1])
