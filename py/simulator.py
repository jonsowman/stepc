# STEPC
# A fast and extensible linear system and controller simulator
# for evaluation of predictive controllers
#
# Jon Sowman 2014 <j.sowman@soton.ac.uk>

import numpy as np


class Simulator:
    timestep = 0
    endtime = 0

    def simulate(self, sys, x0, controller):
        """
        Given a system, an initial system state and a controller, simulate the
        system and return the final system state.

        Examples:
        Simulator::simulate(sys, x0, pidcontroller)
        """

        # Check that parameters are OK
        if(self.timestep == 0):
            raise ValueError("Simulation timestep has not been set")

        if(self.endtime == 0):
            raise ValueError("Simulator end time has not been set")

        if(self.timestep > self.endtime):
            raise ValueError("Simulator end time must be larger than \
                            the timestep")

        # Let the state vector be the initial state
        x = x0

        # Step through time
        i = 0
        while(i < self.endtime):
            # Get the system output
            y = np.dot(sys.C, x)

            # What does the controller want to do
            u = controller.controlmove(y)

            # xdot = Ax + Bu
            xdot = np.dot(sys.A, x)
            xdot += np.dot(sys.B, u)

            # Forwards Euler
            x = x + self.timestep * xdot

            # Increment timestep
            i += self.timestep

        return x
