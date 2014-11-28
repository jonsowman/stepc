/**
 * Fast linear system simulator using Boost
 * and uBLAS.
 *
 * @file controller.hpp
 * @author Jon Sowman, University of Southampton <j.sowman@soton.ac.uk>
 * @copyright Jon Sowman 2014, All Rights Reserved
 * @addtogroup stepc
 */

#pragma once

// Boost
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "linsystem.hpp"

namespace blas = boost::numeric::ublas;

class Controller
{
    public:
        Controller(LinSystem sys)
            : _sys(sys)
        {}

        blas::vector<double> controlmove(const blas::vector<double> y);

    private:
        const LinSystem _sys;
};
