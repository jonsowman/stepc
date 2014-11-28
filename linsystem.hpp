/**
 * Fast linear system simulator using Boost
 * and uBLAS.
 *
 * @file linsystem.hpp
 * @author Jon Sowman, University of Southampton <j.sowman@soton.ac.uk>
 * @copyright Jon Sowman 2014, All Rights Reserved
 * @addtogroup stepc
 */

#pragma once

namespace blas = boost::numeric::ublas;

class LinSystem
{
    public:
        /**
         * Create a linear system
         * @param order The number of system states
         * @param numinputs The number of inputs to the linear system
         * @param numoutputs The number of outputs of the linear system
         */
        LinSystem(const int order, const int numinputs, const int numoutputs)
            : _order(order), _numinputs(numinputs), _numoutputs(numoutputs)
        {
            // Check that the defined system topology is acceptable
            // Having no inputs (an autonomous system) is OK
            if(order == 0 || numoutputs == 0)
            {
                std::cout << "System must have at least 1 state and 1 output"
                    << std::endl;
                return;
            }
            this->A.resize(order, order);
            this->B.resize(order, numinputs);
            this->C.resize(numoutputs, order);
        }

        // Documentation is in the .cpp file
        int getOrder(void);
        int getNumInputs(void);
        int getNumOutputs(void);

        // Define the system matrices
        // Assume D=0 for systems simulated by stepc
        blas::matrix<double> A;
        blas::matrix<double> B;
        blas::matrix<double> C;

    private:
        const int _order, _numinputs, _numoutputs;
};