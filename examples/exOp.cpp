#include <iostream>
#include <unistd.h>
#include <LatSim/Lattice.hpp>

using namespace std;
using namespace LatSim;

int main(int argc, char *argv[])
{
    // Argument parsing
    Coord<4> l, p;
    int      n;

    if (argc != 7)
    {
        cerr << "usage: " << argv[0] << " <L> <p0> <p1> <p2> <p3> <nOp>" << endl;

        return EXIT_FAILURE;
    }
    l[0] = strTo<unsigned int>(argv[1]);
    l[1] = l[0];
    l[2] = l[0];
    l[3] = l[0];
    p[0] = strTo<unsigned int>(argv[2]);
    p[1] = strTo<unsigned int>(argv[3]);
    p[2] = strTo<unsigned int>(argv[4]);
    p[3] = strTo<unsigned int>(argv[5]);
    n    = strTo<int>(argv[6]);

    // Create layout
    Layout<4> layout(l, p);

    registerGlobalLayout(layout);

    // LatSim implementation
    double             time, nOp;
    Lattice<double, 4> f1,f2,f3,f4;

    nOp  = 0.;
    for (unsigned long i = 0; i < layout.getLocalVolume(); ++i)
    {
        f1[i] = 2*static_cast<double>(argc)+i;
        f2[i] = 3*static_cast<double>(argc)+i;
        f3[i] = 1*static_cast<double>(argc)+i;
        f4[i] = sin(i);
    }
    sleep(1);
    time = layout.time();
    for (int j = 0; j < n; ++j)
    {
        f2   = f1+(f2-f3)*f4*f4*f4+f1-f3;
        nOp += 7.*static_cast<double>(layout.getVolume());
    }
    time = layout.time() - time;
    sleep(1);
    if (layout.getRank() == 0)
    {
        cout << time << " s, " << nOp << " flop, " << nOp/time << " flop/s" << endl;
//        for (unsigned long i = 0; i < layout.getLocalVolume()/4; ++i)
//        {
//            cout << "f2[" << i << "]= " << f2[i] << endl;
//        }
    }

    // naive implementation
    nOp  = 0.;
    for (unsigned long i = 0; i < layout.getLocalVolume(); ++i)
    {
        f1[i] = 2*static_cast<double>(argc)+i;
        f2[i] = 3*static_cast<double>(argc)+i;
        f3[i] = 1*static_cast<double>(argc)+i;
        f4[i] = sin(i);
    }
    sleep(1);
    time = layout.time();
    for (int j = 0; j < n; ++j)
    {
        for (unsigned long i = 0; i < layout.getLocalVolume(); ++i)
        {
            f2[i] = f1[i]+(f2[i]-f3[i])*f4[i]*f4[i]*f4[i]+f1[i]-f3[i];
        }
        nOp += 7.*static_cast<double>(layout.getVolume());
    }
    time = layout.time() - time;
    sleep(1);
    if (layout.getRank() == 0)
    {
        cout << time << " s, " << nOp << " flop, " << nOp/time << " flop/s" << endl;
//        for (unsigned long i = 0; i < layout.getLocalVolume()/4; ++i)
//        {
//            cout << "f2[" << i << "]= " << f2[i] << endl;
//        }
    }

    return EXIT_SUCCESS;
}
