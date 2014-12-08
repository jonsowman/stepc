# STEPC
# A fast and extensible linear system and controller simulator
# for evaluation of predictive controllers
#
# Jon Sowman 2014 <j.sowman@soton.ac.uk>

import numpy as np
import subprocess

from stepc.linsystem import LinSystem
from stepc.simulator import Simulator
from stepc.controller import PIDController
from stepc.solver import ForwardsEulerSolver

git_version = subprocess.check_output(["git", "describe", "--dirty",
                                      "--always"])
print "+++++ STEPC (Version %s) +++++" % git_version.strip()

# Create the system
sys = LinSystem(2, 1, 1)

# Mass on a spring
# Mass
m = .5
# Damping const
b = .2
# Spring const
k = 1

sys.A[0, 0] = 0
sys.A[0, 1] = 1
sys.A[1, 0] = -k/m
sys.A[1, 1] = -b/m

sys.B[0, 0] = 0
sys.B[1, 0] = 1/m

sys.C[0, 0] = 1
sys.C[0, 1] = 0

# Now a simple controller
pid = PIDController()
pid.set_kp(2)
# Target state is [0,0]
pid.set_target(np.array([0]))

# Make a simulator and set parms
solver = ForwardsEulerSolver()
sim = Simulator()
sim.set_solver(solver)
sim.set_timestep(0.001)
sim.set_endtime(30)

# Initial condition
x0 = np.array([[-1], [0]])

# Go for it
t_all, x_all, u_all = sim.simulate(sys, x0, pid)
x = x_all[:, -1]

# Results
print "Final state is [%.6f, %.5f]" % (x[0], x[1])
