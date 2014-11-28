/**
 * Fast linear system simulator using Boost
 * and uBLAS.
 *
 * @file linsystem.cpp
 * @author Jon Sowman, University of Southampton <j.sowman@soton.ac.uk>
 * @copyright Jon Sowman 2014, All Rights Reserved
 * @addtogroup stepc
 */

// Boost
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

// Project local
#include "linsystem.hpp"

namespace blas = boost::numeric::ublas;

/**
 * Get order of this system
 * @returns The integer order of this system
 */
int LinSystem::getOrder(void)
{
    return _order;
}

/**
 * Get the number of inputs this system takes
 * @returns The number of inputs
 */
int LinSystem::getNumInputs(void)
{
    return _numinputs;
}

/**
 * Get the number of outputs this system has
 * @returns The number of system outputs
 */
int LinSystem::getNumOutputs(void)
{
    return _numoutputs;
}
