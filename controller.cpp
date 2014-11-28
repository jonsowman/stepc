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
 * This method returns a null control move, since all controllers should be
 * derived classes
 * @warning You should not be using this method!
 * @param y A LinSystem output vector
 * @return A null vector.
 */
blas::vector<double> Controller::controlMove(const blas::vector<double> y)
{
    return y * 0;
}

/**
 * Given a system output y, compute a control move.
 * @param y The current system output
 * @returns A control vector u that will be applied to the system at the next
 * timestep by the Simulator
 */
blas::vector<double> PIDController::controlMove(const blas::vector<double> y)
{
    blas::vector<double> u(1);
    blas::vector<double> e(1);

    // Find error signal
    e = this->_target - y;

    // Add error to the integrator
    this->integrator += e;

    // Calculate the control action in PID fashion
    u = (this->_kp * e);
    u += (this->integrator * this->_ki);
    u += (this->_kd * (e - this->old_error));

    // Store old error
    this->old_error = e;

    // Return the input vector
    return u;
}

void PIDController::setTarget(const blas::vector<double> target)
{
    _target = target;
}

void PIDController::setPropGain(const double kp)
{
    _kp = kp;
}

void PIDController::setIntGain(const double ki)
{
    _ki = ki;
}

void PIDController::setDiffGain(const double kd)
{
    _kd = kd;
}
