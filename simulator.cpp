/**
 * Fast linear system simulator using Boost
 * and uBLAS.
 *
 * @file simulator.cpp
 * @author Jon Sowman, University of Southampton <j.sowman@soton.ac.uk>
 * @copyright Jon Sowman 2014, All Rights Reserved
 * @addtogroup stepc
 */

#include <iostream>
#include <vector>

// Boost
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "simulator.hpp"

namespace blas = boost::numeric::ublas;

/**
 * Implementation of the linear system simulator
 * @param sys A LinSystem to be simulated
 * @param x The initial state vector for the system to be simulated
 * @returns The final state
 */
blas::vector<double> Simulator::simulate(LinSystem sys, blas::vector<double> x)
{
    // Sanity check simulation parameters
    if(_ts == 0.0 || _endtime == 0.0 || _ts > _endtime)
    {
        std::cout << "Sim params don't make sense!" << std::endl;
        return x;
    }

    // Step in time
    for(double i = 0.0; i<=_endtime; i+=_ts)
    {
        // FIXME u should be found byt he controller
        blas::vector<double> u(1);
        u(0) = 0.1;
        blas::vector<double> xdot = prod(sys.A, x) + prod(sys.B, u);

        // Integrate using forwards Euler (could do better here)
        x = x + _ts * xdot;
    }
    return x;
}
