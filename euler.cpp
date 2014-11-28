/**
 * Fast linear system simulator using Boost
 * and uBLAS.
 *
 * @file euler.cpp
 * @author Jon Sowman, University of Southampton <j.sowman@soton.ac.uk>
 * @copyright Jon Sowman 2014, All Rights Reserved
 * @addtogroup stepc
 */

#include <iostream>
#include <vector>

// Boost
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

// Project local
#include "simulator.hpp"
#include "linsystem.hpp"

namespace blas = boost::numeric::ublas;

int main(void)
{
    std::cout << "------ STEPC ------" << std::endl;
    blas::vector<double> x(2);

    x(0) = 1.0;
    x(1) = 0.0;

    // Create linear system of order 2 and 1 input and 1 output
    LinSystem sys(2, 1, 1);

    // Define A (2x2)
    sys.A(0,0) = 0;
    sys.A(0,1) = 1;
    sys.A(1,0) = -1;
    sys.A(1,1) = 0;

    // Define B (2x1)
    sys.B(0,0) = 0;
    sys.B(1,0) = -1;

    // Define C (1x2)
    sys.C(0,0) = 1;
    sys.C(0,1) = 0;

    // Create a fixed timestep simulator
    Simulator Sim;
    Sim.setTimestep(0.001);
    Sim.setEndTime(10);

    // Simulate 'sys' from our initial state 'x'
    x = Sim.simulate(sys, x);

    // Final
    std::cout << "Final state is: [" << x(0) 
        << ", " << x(1) << "] " << std::endl;

    return 0;
}

