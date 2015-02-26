/*
 * Global.hpp, part of LatSim
 *
 * Copyright (C) 2014-2015 Antonin Portelli
 *
 * LatSim is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * LatSim is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with LatSim.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef LatSim_Global_hpp_
#define	LatSim_Global_hpp_

#include <LatCore/LatCore.hpp>
#include <mpi.h>

#ifdef LATSIM_FORCE_INLINE
#define strong_inline __attribute__((always_inline)) inline
#else
#define strong_inline inline
#endif

#ifdef LATSIM_FLATTEN
#define flatten __attribute__((flatten))
#else
#define flatten
#endif

#define BEGIN_LATSIM_NAMESPACE \
namespace LatSim {\
using namespace LatCore;
#define END_LATSIM_NAMESPACE }

BEGIN_LATSIM_NAMESPACE

// lattice coordinate type /////////////////////////////////////////////////////
template <unsigned int D>
using Coord = std::array<unsigned int, D>;

// lattice base type ///////////////////////////////////////////////////////////
class LatticeObj {};

// Error handling //////////////////////////////////////////////////////////////
#define LATSIM_SRC_LOC LatSim::strFrom(__FUNCTION__) + " at " +\
                LatSim::strFrom(__FILE__) + ":" + LatSim::strFrom(__LINE__)
#define locGlobalError(msg) LatSim::globalError(msg, LATSIM_SRC_LOC)

void inline globalError(const std::string msg, const std::string loc = "")
{
    int rank, isInit;

    MPI_Initialized(&isInit);
    if (!isInit)
    {
        MPI_Init(0, nullptr);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        std::cerr << "GLOBAL ERROR: " << msg;
        if (!loc.empty())
        {
            std::cerr << " (" << loc << ")" << std::endl;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
    {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
}

// Indexing helpers ////////////////////////////////////////////////////////////
template <size_t D>
inline unsigned int coordToRowMajor(const Coord<D> &x, const Coord<D> &dim)
{
    unsigned int ind;

    ind = dim[1]*x[0];
    for (unsigned int d = 1; d < D - 1; ++d)
    {
        ind = dim[d+1]*(x[d] + ind);
    }
    ind += x[D-1];

    return ind;
}

template <size_t D>
inline Coord<D> rowMajorToCoord(const unsigned int ind, const Coord<D> &dim)
{
    Coord<D>     x;
    unsigned int j, dimprod, ud;

    j       = ind;
    dimprod = 1;
    for (int d = D - 1; d >= 0; --d)
    {
        ud       = static_cast<unsigned int>(d);
        x[ud]    = (j/dimprod)%dim[ud];
        j       -= dimprod*x[ud];
        dimprod *= dim[ud];
    }

    return x;
}

END_LATSIM_NAMESPACE

#endif // LatSim_Global_hpp_
