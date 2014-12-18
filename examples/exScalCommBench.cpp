#include <iostream>
#include <sstream>
#include <unistd.h>
#include <LatSim/Lattice.hpp>
#include <LatSim/LatOp.hpp>

using namespace std;
using namespace LatSim;

const unsigned int sleepTime = 1;

class CommBench: public Logger
{
public:
    CommBench(const int n): Logger("CommBench"), n_(n){};
    ~CommBench(void) = default;
    void run(void);
private:
    int n_;
};

void CommBench::run(void)
{
    Layout<4>         &layout = *dynamic_cast<Layout<4> *>(globalLayout);
    Lattice<float, 4> x, y, z;
    stringstream      buf;
    double            time, nOp, nComm;

    for (unsigned int i = 0; i < layout.getLocalVolume(); ++i)
    {
        y(i) = 3.f/layout.getDim(0)+i;
    }

    masterLog("* site-wise subtraction...");
    sleep(sleepTime);
    time = layout.time();
    for (int j = 0; j < n_; ++j)
    {
        z = x - y;
    }
    time = layout.time() - time;
    sleep(sleepTime);
    nOp = n_*static_cast<float>(layout.getVolume());
    buf << time << " s, " << nOp << " flop, " << nOp/time << " flop/s 0 MB/s";
    masterLog(buf.str());

    masterLog("* Laplacian...");
    sleep(sleepTime);
    time = layout.time();
    for (int j = 0; j < n_; ++j)
    {
        x = laplacian(y);
    }
    time = layout.time() - time;
    sleep(sleepTime);
    nOp   = n_*9.*static_cast<float>(layout.getVolume());
    nComm = 0.;
    for (unsigned int d = 0; d < 8; ++d)
    {
        if (layout.needComm(d))
        {
            nComm += n_*static_cast<float>(sizeof(float)
                                          *layout.getLocalSurface(d));
        }
    }
    buf.str("");
    buf << time << " s, " << nOp << " flop, " << nOp/time << " flop/s";
    buf << " " << nComm/time/1024./1024. << " MB/s";
    masterLog(buf.str());
    masterLog("* forward finite differences...");
    for (unsigned int d = 0; d < 8; ++d)
    {

        sleep(sleepTime);
        time = layout.time();
        for (int j = 0; j < n_; ++j)
        {
            x = forwardDiff(y, d);
        }
        time = layout.time() - time;
        sleep(sleepTime);
        nOp = n_*static_cast<float>(layout.getVolume());
        if (layout.needComm(d))
        {
            nComm = n_*static_cast<float>(sizeof(float)
                                         *layout.getLocalSurface(d));
        }
        else
        {
            nComm = 0.;
        }
        buf.str("");
        buf << time << " s, " << nOp << " flop, " << nOp/time << " flop/s";
        buf << " " << nComm/time/1024./1024. << " MB/s";
        masterLog(buf.str());
    }
}

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

    // Benchmark
    CommBench bench(n);

    bench.run();

    return EXIT_SUCCESS;
}
