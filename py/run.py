# STEPC
# A fast and extensible linear system and controller simulator
# for evaluation of predictive controllers
#
# Jon Sowman 2014 <j.sowman@soton.ac.uk>

import numpy as np
import subprocess

from linsystem import LinSystem
from simulator import Simulator
from controller import PIDController

git_version = subprocess.check_output(["git", "describe", "--dirty",
                                      "--always"])
print "+++++ STEPC (Version %s) +++++" % git_version.strip()

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

# Now a simple controller
pid = PIDController()
pid.set_kp(0.1)
# Target state is [0,0]
pid.set_target(np.matrix(0))

# Make a simulator and set parms
sim = Simulator()
sim.set_timestep(0.001)
sim.set_endtime(10)

# Initial condition
x0 = np.matrix([[-1], [0]])

# Go for it
t_all, x_all, u_all = sim.simulate(sys, x0, pid)
x = x_all[:, -1]

# Results
print "Final state is [%.6f, %.5f]" % (x[0], x[1])
