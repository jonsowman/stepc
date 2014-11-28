/**
 * Fast linear system simulator using Boost
 * and uBLAS.
 *
 * @file simulator.hpp
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

class Simulator
{
    public:

    void setTimestep(double Ts)
    {
        if(Ts != 0.0)
            _ts = Ts;
    }

    void setEndTime(double EndTime)
    {
        if(EndTime != 0.0)
            _endtime = EndTime;
    }

    blas::vector<double> simulate(LinSystem sys, blas::vector<double> x);

    private:
        double _ts;
        double _endtime;
};

