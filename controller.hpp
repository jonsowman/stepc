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
            this->_target.resize(sys->getNumOutputs());
            this->integrator.resize(sys->getNumOutputs());
            this->old_error.resize(sys->getNumOutputs());
        }

        // Destroy a PIDController
        ~PIDController() {}

        virtual blas::vector<double> controlMove(const blas::vector<double> y);

        void setTarget(const blas::vector<double> target);
        void setPropGain(const double kp);
        void setIntGain(const double ki);
        void setDiffGain(const double kd);

    private:
        /// A reference vector for the inputs
        blas::vector<double> _target;
        
        /// K_p gains
        double _kp;

        /// K_i gains
        double _ki;

        /// K_d gains
        double _kd;

        /// Integrator for K_i gain
        blas::vector<double> integrator;

        /// Retain one past error for the K_d term
        blas::vector<double> old_error;
};
