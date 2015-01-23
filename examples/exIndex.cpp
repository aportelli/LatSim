#include <iostream>
#include <unistd.h>

#include <LatSim/Lattice.hpp>

#define L 4
#define D 2

using namespace std;
using namespace LatSim;

int main(void)
{
    // Create layout
    Coord<4>     l{{L, L, L, L}}, p{{1, 1, 1, 1}}, x, xn;

    p[D]  = 2;
    l[D] *= 2;

    Layout<4>    layout(l, p);

    registerGlobalLayout(layout);

    // Dump nearest neighbour table
    Coord<4> ld = layout.getLocalDim();

    if (layout.getMyRank() == 0)
    {
        cout << "-- NEAREST NEIGHBOUR TABLE" << endl;
        for (unsigned int i = 0; i < layout.getLocalVolume(); ++i)
        {
            x = rowMajorToCoord<4>(i, ld);

            cout << "[";
            for (unsigned int d = 0; d < 4; ++d)
            {
                cout << x[d] << ((d == 3) ? "" : ", ");
            }
            cout << "] -> " << i << endl;

            for (unsigned int d = 0; d < 8; ++d)
            {
                cout << "NN (d= " << d << ") -> ";
                cout << layout.getNearNeigh(i, d) << endl;
            }
        }
    }

    // Check communications are mapped correctly
    Lattice<unsigned int, 4> lat;

    for (unsigned int i = 0; i < layout.getLocalVolume(); ++i)
    {
        lat(i) = i;
    }
    for (unsigned int d = 0; d < 8; ++d)
    {
        lat.gather(d);
        if (layout.getMyRank() == 0)
        {
            cout << "-- NEIGHBOUR IN DIRECTION " << d << endl;
            for (unsigned int i = 0; i < layout.getLocalVolume(); ++i)
            {
                x  = rowMajorToCoord(i, ld);
                xn = rowMajorToCoord(lat(i, d), ld);
                cout << "lat([";
                for (unsigned int d = 0; d < 4; ++d)
                {
                    cout << x[d] << ((d == 3) ? "" : ", ");
                }
                cout << "](" << i << "), " << d << ")= ";
                cout << "[";
                for (unsigned int d = 0; d < 4; ++d)
                {
                    cout << xn[d] << ((d == 3) ? "" : ", ");
                }
                cout << "](" << lat(i, d) << ")" << endl;
            }
        }
    }

    return EXIT_SUCCESS;
}
