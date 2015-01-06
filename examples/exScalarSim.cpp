#include <iostream>
#include <LatSim/Lattice.hpp>
#include <LatSim/VerletIntegrator.hpp>
#include <LatSim/OmelyanIntegrator.hpp>

using namespace std;
using namespace LatSim;

typedef Lattice<double, 4> LScalar4;

int main(void)
{
    VerletIntegrator<LScalar4, LScalar4> integrator(nullptr, nullptr);
    LScalar4 phi, pi;

    integrator.evolve(phi, pi, 0.1, 10);
    
    return 0;
}
