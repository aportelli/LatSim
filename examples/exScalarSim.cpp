#include <iostream>
#include <LatSim/Lattice.hpp>
#include <LatSim/Hmc.hpp>
#include <LatSim/VerletIntegrator.hpp>
#include <LatSim/OmelyanIntegrator.hpp>
#include <LatSim/Phi4.hpp>

using namespace std;
using namespace LatSim;
using namespace Phi4;

constexpr unsigned int D = 4;

typedef Lattice<double, D> LScalar;

int main(int argc, char *argv[])
{
    Coord<D>     l, p;
    ActionParam  param;
    unsigned int nTraj, nStep;

    if (argc != 10)
    {
        cerr << "usage: " << argv[0];
        cerr << " <L> <p0> <p1> <p2> <p3> <lambda> <kappa> <nTraj> <nStep>";
        cerr << endl;

        return EXIT_FAILURE;
    }
    l[0]         = strTo<unsigned int>(argv[1]);
    l[1]         = l[0];
    l[2]         = l[0];
    l[3]         = l[0];
    p[0]         = strTo<unsigned int>(argv[2]);
    p[1]         = strTo<unsigned int>(argv[3]);
    p[2]         = strTo<unsigned int>(argv[4]);
    p[3]         = strTo<unsigned int>(argv[5]);
    param.lambda = strTo<double>(argv[6]);
    param.kappa  = strTo<double>(argv[7]);
    nTraj        = strTo<unsigned int>(argv[8]);
    nStep        = strTo<unsigned int>(argv[9]);

    Layout<D> layout(l, p);

    registerGlobalLayout(layout);

    MoveMom<D>                          moveMom(param);
    Hamiltonian<D>                      h(param);
    OmelyanIntegrator<LScalar, LScalar> integrator(moveField<D>, moveMom);
    Hmc<LScalar, LScalar>               hmc(h, momGen<D>, integrator, nStep);
    LScalar                             phi, pi;

    phi.fill(0.0);
    for (unsigned int t = 0; t < nTraj; ++t)
    {
        double phiAv, locSum = 0.;

        if(hmc.update(phi, pi))
        {
            FOR_SITE(phi, i)
            {
                locSum += phi(i);
            }
            MPI_Reduce(&locSum, &phiAv, 1, MPI_DOUBLE, MPI_SUM, 0,
                       layout.getCommGrid());
            if (layout.getRank() == 0)
            {
                cout << "<phi>= " << phiAv/static_cast<double>(layout.getVolume());
                cout << endl;
            }
        }
    }
    hmc.printStat();
    
    return 0;
}
