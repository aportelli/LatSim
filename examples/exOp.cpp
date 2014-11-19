#include <iostream>
#include <LatSim/Lattice.hpp>

using namespace std;
using namespace LatSim;

int main(void)
{
    Coord<4>  l{{32,32,32,32}}, p{{2,2,2,1}};
    Layout<4> layout(l, p);

    registerGlobalLayout(layout);

    Lattice<double, 4> phi;

    phi.gather(4);
    
    return EXIT_SUCCESS;
}
