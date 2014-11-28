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
        Controller(LinSystem *sys)
            : _sys(sys)
        {}

        virtual ~Controller() {}

        // Different controllers implement this differently
        virtual blas::vector<double> controlMove(const blas::vector<double> y);

    private:
        const LinSystem *_sys;
};

class PIDController: public Controller
{
    public:
        // Create a PIDController and call the Controller constructor
        PIDController(LinSystem *sys)
        : Controller(sys)
        {
            // Size the target vector
            _target.resize(sys->getNumOutputs());
        }

        // Destroy a PIDController
        ~PIDController() {}

        virtual blas::vector<double> controlMove(const blas::vector<double> y);

        void setTarget(const blas::vector<double> target);
        void setPropGain(const blas::vector<double> kp);

    private:
        /// A reference vector for the inputs
        blas::vector<double> _target;
        
        /// A vector of K_p gains
        blas::vector<double> _kp;
};
