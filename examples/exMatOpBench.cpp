#include <iostream>
#include <unistd.h>

#include <LatSim/Lattice.hpp>

#define N 16

using namespace std;
using namespace LatSim;

const unsigned int sleepTime = 1;

#define func(out, x, y, z) out = x*(z*(y-x)+y*z*z*(x-y)+z+x)
#define mulFlop (N*N*(2.*N-1.))
#define cwFlop (N*N)
#define funcFlop (5.*mulFlop+5.*cwFlop)

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
    double                  time, nOp;
    Lattice<SFMat<N, N>, 4> x, y, z;

    for (unsigned int i = 0; i < layout.getLocalVolume(); ++i)
    {
        x[i].fill(2*static_cast<float>(argc)+i);
        y[i].fill(3*static_cast<float>(argc)+i);
        z[i].fill(sin(static_cast<float>(i)*argc));
    }
    sleep(sleepTime);
    time = layout.time();
    for (int j = 0; j < n; ++j)
    {
        func(z, x, y, z);
    }
    time = layout.time() - time;
    sleep(sleepTime);
    nOp = n*funcFlop*static_cast<float>(layout.getVolume());
    if (layout.getRank() == 0)
    {
        cout << time << " s, " << nOp << " flop, " << nOp/time << " flop/s" << endl;
    }

    // naive implementation with lattices
    for (unsigned int i = 0; i < layout.getLocalVolume(); ++i)
    {
        x[i].fill(2*static_cast<float>(argc)+i);
        y[i].fill(3*static_cast<float>(argc)+i);
        z[i].fill(sin(static_cast<float>(i)*argc));
    }
    sleep(sleepTime);
    time = layout.time();
    for (int j = 0; j < n; ++j)
    {
        for (unsigned int i = 0; i < layout.getLocalVolume(); ++i)
        {
            func(z[i], x[i], y[i], z[i]);
        }
    }
    time = layout.time() - time;
    sleep(sleepTime);
    if (layout.getRank() == 0)
    {
        cout << time << " s, " << nOp << " flop, " << nOp/time << " flop/s" << endl;
    }

    // naive implementation with arrays
    SFMat<N, N> *a = new SFMat<N, N>[layout.getLocalVolume()];
    SFMat<N, N> *b = new SFMat<N, N>[layout.getLocalVolume()];
    SFMat<N, N> *c = new SFMat<N, N>[layout.getLocalVolume()];

    for (unsigned int i = 0; i < layout.getLocalVolume(); ++i)
    {
        a[i].fill(2*static_cast<float>(argc)+i);
        b[i].fill(3*static_cast<float>(argc)+i);
        c[i].fill(sin(static_cast<float>(i)*argc));
    }
    sleep(sleepTime);
    time = layout.time();
    for (int j = 0; j < n; ++j)
    {
        for (unsigned int i = 0; i < layout.getLocalVolume(); ++i)
        {
            func(c[i], a[i], b[i], c[i]);
        }
    }
    time = layout.time() - time;
    sleep(sleepTime);
    if (layout.getRank() == 0)
    {
        cout << time << " s, " << nOp << " flop, " << nOp/time << " flop/s" << endl;
    }

    delete [] a;
    delete [] b;
    delete [] c;

    return EXIT_SUCCESS;
}
