/**
 * Fast linear system simulator using Boost
 * and uBLAS.
 *
 * Jon Sowman 2014 <jon@jonsowman.com>
 *
 * All Rights Reserved
 */

#include <iostream>
#include <vector>

// Boost
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace blas = boost::numeric::ublas;

// Class prototypes
class Simulator;
class LinSystem;

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

class LinSystem
{
    public:
        LinSystem(const int order)
        {
            A.resize(order, order);
        }

        blas::matrix<double> A;
};

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
        // xdot = A*x
        blas::vector<double> xdot = prod(sys.A, x);

        // Integrate using forwards Euler (could do better here)
        x = x + _ts * xdot;
    }
    return x;
}

int main(void)
{
    std::cout << "------ STEPC ------" << std::endl;
    blas::vector<double> x(2);

    x(0) = 1.0;
    x(1) = 0.0;

    // Create linear system of order 2
    LinSystem sys(2);
    sys.A(0,0) = 0;
    sys.A(0,1) = 1;
    sys.A(1,0) = -1;
    sys.A(1,1) = 0;

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

