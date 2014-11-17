/*
 * Layout.hpp, part of LatSim
 *
 * Copyright (C) 2014 Antonin Portelli
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

#ifndef LatSim_Layout_hpp_
#define LatSim_Layout_hpp_

#include <LatSim/Global.hpp>
#include <LatSim/Logger.hpp>

BEGIN_NAMESPACE

/******************************************************************************
 *                                 Layout                                     *
 ******************************************************************************/
class LayoutObject
{
public:
    LayoutObject(void)          = default;
    virtual ~LayoutObject(void) = default;
};

// layout class
template <unsigned long D>
class Layout: public LayoutObject, public Logger
{
public:
    // constructor
    Layout(const Coord<D> &l, const Coord<D> &p);
    // destructor
    virtual ~Layout(void);
    // access
    unsigned int  getDim(const unsigned int d) const;
    unsigned int  getLocalDim(const unsigned int d) const;
    unsigned long getLocalSurface(const unsigned int d) const;
    unsigned long getVolume(void) const;
    unsigned long getLocalVolume(void) const;
    // check
    static void checkDir(const unsigned int d);
    // get absolute direction
    static unsigned int absDir(const unsigned int d);
    // get opposite direction
    static unsigned int oppDir(const unsigned int d);
    // get neighbor coordinate in the grid
    int neighborCoord(const unsigned int d);
private:
    int                          rank_;
    unsigned long                volume_, locVolume_;
    std::array<unsigned long, D> locSurface_;
    std::array<int, D>           p_, coord_;
    Coord<D>                     dim_, locDim_;
    MPI_Comm                     commGrid_;
    std::array<MPI_Comm, D>      commDirGrid_;
};

// global layout
extern LayoutObject *globalLayout;

template <unsigned long D>
void registerGlobalLayout(Layout<D> &layout)
{
    globalLayout = &layout;
}

/******************************************************************************
 *                           Layout implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <unsigned long D>
Layout<D>::Layout(const Coord<D> &dim, const Coord<D> &p)
: Logger("Layout")
, dim_(dim)
{
    int         size, isInit, nproc = 1, isPeriodic[D], isActive[D];
    std::string buf;

    // MPI initialization
    MPI_Initialized(&isInit);
    if (isInit)
    {
        locGlobalError("MPI already initialized (it was manually initialized"
                       " or several Layout objects were created)");
    }
    MPI_Init(0, nullptr);

    // dimension check
    for (unsigned int d = 0; d < D; ++d)
    {
        p_[d]   = static_cast<int>(p[d]);
    }
    buf = "";
    for (unsigned int d = 0; d < D; ++d)
    {
        buf += strFrom(dim_[d]) + ((d != D - 1) ? "x" : "");
    }
    masterLog("lattice size   : " + buf);
    buf = "";
    for (unsigned int d = 0; d < D; ++d)
    {
        buf += strFrom(p_[d]) + ((d != D - 1) ? "x" : "");
    }
    masterLog("MPI partition  : " + buf);
    for (unsigned int d = 0; d < D; ++d)
    {
        if (dim_[d] % p_[d] != 0)
        {
            locGlobalError("lattice size (" + strFrom(dim_[d]) +
                           ") is not divisible by the number of processes (" +
                           strFrom(p_[d]) + ") in direction " + strFrom(d));
        }
    }
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    for (unsigned int d = 0; d < D; ++d)
    {
        nproc *= p_[d];
    }
    if (nproc != size)
    {
        locGlobalError("layout number of processes (" + strFrom(nproc) +
                       ") does not match the real one (" + strFrom(size) + ")");
    }

    // compute volumes and surfaces
    volume_    = 1;
    locVolume_ = 1;
    for (unsigned int d = 0; d < D; ++d)
    {
        locDim_[d]      = dim_[d]/p_[d];
        volume_        *= dim_[d];
        locVolume_     *= locDim_[d];
        locSurface_[d]  = 1;
        for (unsigned int i = 0; i < D; ++i)
        {
            locSurface_[d] *= (i != d) ? locDim_[d] : 1;
        }
    }
    masterLog("volume         : " + strFrom(volume_));
    masterLog("local volume   : " + strFrom(locVolume_));
    buf = "";
    for (unsigned int d = 0; d < D; ++d)
    {
        buf += strFrom(locSurface_[d]) + " (d= " + strFrom(d) + ")";
        buf += (d == D - 1) ? "" : ", ";
    }
    masterLog("local surfaces : " + buf);

    // create Cartesian topology
    for (unsigned int d = 0; d < D; ++d)
    {
        isPeriodic[d] = 1;
    }
    MPI_Cart_create(MPI_COMM_WORLD, D, p_.data(), isPeriodic, 1, &commGrid_);
    for (unsigned int d = 0; d < D; ++d)
    {
        for (unsigned int i = 0; i < D; ++i)
        {
            isActive[i] = (i == d) ? 1 : 0;
        }
        MPI_Cart_sub(commGrid_, isActive, &commDirGrid_[d]);
    }
    MPI_Comm_rank(commGrid_, &rank_);
    MPI_Cart_coords(commGrid_, rank_, D, coord_.data());

    // everything is ok
    masterLog("communication grid ready");
}

// destructor //////////////////////////////////////////////////////////////////
template <unsigned long D>
Layout<D>::~Layout(void)
{
    MPI_Finalize();
    masterLog("communication grid terminated");
}

// access //////////////////////////////////////////////////////////////////////
template <unsigned long D>
unsigned int Layout<D>::getDim(const unsigned int d) const
{
    checkDir(d);

    return dim_[absDir(d)];
}

template <unsigned long D>
unsigned int Layout<D>::getLocalDim(const unsigned int d) const
{
    checkDir(d);

    return locDim_[absDir(d)];
}

template <unsigned long D>
unsigned long Layout<D>::getLocalSurface(const unsigned int d) const
{
    checkDir(d);

    return locSurface_[absDir(d)];
}

template <unsigned long D>
unsigned long Layout<D>::getVolume(void) const
{
    return volume_;
}

template <unsigned long D>
unsigned long Layout<D>::getLocalVolume(void) const
{
    return locVolume_;
}

// check ///////////////////////////////////////////////////////////////////////
template <unsigned long D>
void Layout<D>::checkDir(const unsigned int d)
{
    if (d >= 2*D)
    {
        locGlobalError("direction " + strFrom(d) + " invalid " +
                       "(number of dimension is " + strFrom(D) + ")");
    }
}

// get absolute direction //////////////////////////////////////////////////////
template <unsigned long D>
unsigned int Layout<D>::absDir(const unsigned int d)
{
    return d % D;
}

// get opposite direction //////////////////////////////////////////////////////
template <unsigned long D>
unsigned int Layout<D>::oppDir(const unsigned int d)
{
    return (d + D) % (2*D);
}

// get neighbor coordinate in the grid /////////////////////////////////////////
template <unsigned long D>
int Layout<D>::neighborCoord(const unsigned int d)
{
    checkDir(d);

    return (d < D) ? ((coord_[d] + 1)         % p_[d])
                   : ((coord_[d] - 1 + p_[d]) % p_[d]);
}

END_NAMESPACE

#endif // LatSim_Layout_hpp_
