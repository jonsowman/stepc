/**
 * Fast linear system simulator using Boost
 * and uBLAS.
 *
 * Implements a set of controllers that can be attached to a linear system
 * LinSystem to control its output(s) to their desired value.
 *
 * @file controller.cpp
 * @author Jon Sowman, University of Southampton <j.sowman@soton.ac.uk>
 * @copyright Jon Sowman 2014, All Rights Reserved
 * @addtogroup stepc
 */

// Boost
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "controller.hpp"

namespace blas = boost::numeric::ublas;

/**
 * Given a system output y, compute a control move.
 * @param y The current system output
 * @returns A control vector u that will be applied to the system at the next
 * timestep by the Simulator
 */
blas::vector<double> controlmove(const blas::vector<double> y)
{
    blas::vector<double> x;
    x = 0.1 * y;
    return x;
}
