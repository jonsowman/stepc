# STEPC
# A fast and extensible linear system and controller simulator
# for evaluation of predictive controllers
#
# Jon Sowman 2014 <j.sowman@soton.ac.uk>

import numpy as np


class Simulator(object):
    __timestep = 0
    __endtime = 0
    __solver = False

    def simulate(self, sys, x0, controller):
        """
        Given a system, an initial system state and a controller, simulate the
        system and return the final system state.

        Examples:
        Simulator::simulate(sys, x0, pidcontroller)
        """

        # Did someone set a solver
        if not self.__solver:
            raise AssertionError("No solver set in Simulator!")

        # Check that parameters are OK
        if(self.__timestep == 0):
            raise ValueError("Simulation timestep has not been set")

        if(self.__endtime == 0):
            raise ValueError("Simulator end time has not been set")

        if(self.__timestep > self.__endtime):
            raise ValueError("Simulator end time must be larger than \
                            the timestep")

        # Let the state vector be the initial state
        x = x0

        # Store the system state and inputs
        t_store = np.empty(0)
        x_store = np.empty([sys.order, 0])
        u_store = np.empty([sys.numinputs, 0])

        # Step through time
        i = 0
        while(i < self.__endtime):
            # Get the system output
            y = np.dot(sys.C, x)

            # What does the controller want to do
            u = controller.controlmove(y)

            # xdot = Ax + Bu
            xdot = np.dot(sys.A, x)
            xdot += np.dot(sys.B, u)

            # Call whatever solver we're using
            x = self.__solver.step(x, xdot, self.__timestep)

            # Increment timestep
            i += self.__timestep

            # Store the values
            t_store = np.append(t_store, i)
            x_store = np.append(x_store, x, 1)
            u_store = np.append(u_store, u, 1)

        return (t_store, x_store, u_store)

    def set_timestep(self, timestep):
        self.__timestep = timestep

    def set_endtime(self, endtime):
        self.__endtime = endtime

    def set_solver(self, solver):
        self.__solver = solver
