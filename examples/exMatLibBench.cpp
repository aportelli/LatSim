#include <iomanip>
#include <iostream>
#include <unistd.h>

#include <LatSim/Global.hpp>

using namespace std;
using namespace LatSim;

template <unsigned int N>
void matMulBench(const unsigned long nMul)
{
    SFMat<N, N>      m = SFMat<N, N>::Random(), n = SFMat<N, N>::Random(), p;
    constexpr double nOp = N*N*(2.*N - 1.);
    double           time, avTime;
    int              rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    time = MPI_Wtime();
    for (unsigned int i = 0; i < nMul; ++i)
    {
        m = m*n;
    }
    time = MPI_Wtime() - time;
    MPI_Reduce(&time, &avTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        avTime /= size;
        cout << setw(10) << strFrom(N) + "x" + strFrom(N) + " MUL:";
        cout << setw(20) << strFrom(time) + "s";
        cout << setw(20) << strFrom(size*nMul*nOp) + " flop";
        cout << setw(20) << strFrom((size*nMul*nOp)/time) + " flop/s" << endl;
    }
}

template <unsigned int N>
void matAddBench(const unsigned long nAdd)
{
    SFMat<N, N>      m = SFMat<N, N>::Random(), n = SFMat<N, N>::Random(), p;
    constexpr double nOp = N*N;
    double           time, avTime;
    int              rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    time = MPI_Wtime();
    for (unsigned int i = 0; i < nAdd; ++i)
    {
        m = m+n;
    }
    time = MPI_Wtime() - time;
    MPI_Reduce(&time, &avTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        avTime /= size;
        cout << setw(10) << strFrom(N) + "x" + strFrom(N) + " ADD:";
        cout << setw(20) << strFrom(time) + "s";
        cout << setw(20) << strFrom(size*nAdd*nOp) + " flop";
        cout << setw(20) << strFrom((size*nAdd*nOp)/time) + " flop/s" << endl;
    }
}

template <unsigned int N>
void matExprBench(const unsigned long nExpr)
{
    SFMat<N, N>      x = SFMat<N, N>::Random(), y = SFMat<N, N>::Random(),
                     z = SFMat<N, N>::Random();
    constexpr double nOp = 5.*N*N + 5.*N*N*(2.*N - 1.);
    double           time, avTime;
    int              rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    time = MPI_Wtime();
    for (unsigned int i = 0; i < nExpr; ++i)
    {
        z = x*(z*(y-x)+y*z*z*(x-y)+z+x);
    }
    time = MPI_Wtime() - time;
    MPI_Reduce(&time, &avTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        avTime /= size;
        cout << setw(10) << strFrom(N) + "x" + strFrom(N) + " EXPR:";
        cout << setw(20) << strFrom(time) + "s";
        cout << setw(20) << strFrom(size*nExpr*nOp) + " flop";
        cout << setw(20) << strFrom((size*nExpr*nOp)/time) + " flop/s" << endl;
    }
}

int main(int argc __unused, char *argv[])
{
    unsigned int nOp = strTo<unsigned int>(argv[1]);

    MPI_Init(0, nullptr);

    matMulBench<2>(nOp);
    sleep(1);
    matMulBench<3>(nOp);
    sleep(1);
    matMulBench<4>(nOp);
    sleep(1);
    matMulBench<5>(nOp);
    sleep(1);
    matMulBench<6>(nOp);
    sleep(1);
    matMulBench<7>(nOp);
    sleep(1);
    matMulBench<8>(nOp);
    sleep(1);
    matMulBench<9>(nOp);
    sleep(1);

    matAddBench<2>(nOp);
    sleep(1);
    matAddBench<3>(nOp);
    sleep(1);
    matAddBench<4>(nOp);
    sleep(1);
    matAddBench<5>(nOp);
    sleep(1);
    matAddBench<6>(nOp);
    sleep(1);
    matAddBench<7>(nOp);
    sleep(1);
    matAddBench<8>(nOp);
    sleep(1);
    matAddBench<9>(nOp);
    sleep(1);

    matExprBench<2>(nOp);
    sleep(1);
    matExprBench<3>(nOp);
    sleep(1);
    matExprBench<4>(nOp);
    sleep(1);
    matExprBench<5>(nOp);
    sleep(1);
    matExprBench<6>(nOp);
    sleep(1);
    matExprBench<7>(nOp);
    sleep(1);
    matExprBench<8>(nOp);
    sleep(1);
    matExprBench<9>(nOp);

    MPI_Finalize();
}
