# STEPC
# A fast and extensible linear system and controller simulator
# for evaluation of predictive controllers
#
# Jon Sowman 2014 <j.sowman@soton.ac.uk>

import numpy as np

class Simulator:

    timestep = 0
    endtime = 0

    def simulate(self, sys, x0):
        # Check that parameters are OK
        if(self.timestep == 0 or self.endtime == 0):
            print "Simulation timestep or end time are wrong"
            return

        # Let the state vector be the initial state
        x = x0

        # Step through time
        i = 0
        while(i < self.endtime):
            # Get the system output
            y = np.dot(sys.C, x);

            # xdot = Ax + Bu
            xdot = np.dot(sys.A, x)

            # Forwards Euler
            x = x + self.timestep * xdot

            # Increment timestep
            i += self.timestep

        return x
