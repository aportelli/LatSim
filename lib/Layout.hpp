/*
 * Layout.hpp, part of LatSim
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

#ifndef LatSim_Layout_hpp_
#define LatSim_Layout_hpp_

#include <LatSim/Global.hpp>
#include <LatSim/Logger.hpp>
#include <random>

#ifndef RNGTYPE
#define RNGTYPE std::mt19937
#endif

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
template <unsigned int D>
class Layout: public LayoutObject, public Logger
{
public:
    struct PlaneInfo
    {
        unsigned int nBlocks, blockSize, stride;
    };
    typedef RNGTYPE              RngType;
    typedef RngType::result_type SeedType;
public:
    // constructor
    Layout(const Coord<D> &l, const Coord<D> &p);
    // destructor
    virtual ~Layout(void);
    // access
    unsigned int         getRank(void) const;
    const Coord<D> &     getDim(void) const;
    unsigned int         getDim(const unsigned int d) const;
    const Coord<D> &     getLocalDim(void) const;
    unsigned int         getLocalDim(const unsigned int d) const;
    const Coord<D> &     getLocalSurface(void) const;
    unsigned int         getLocalSurface(const unsigned int d) const;
    unsigned int         getVolume(void) const;
    unsigned int         getLocalVolume(void) const;
    unsigned int         getCommBufferSize(void) const;
    unsigned int         getNearNeigh(const unsigned int i,
                                      const unsigned int d) const;
    const unsigned int * getNearNeigh(const unsigned int d) const;
    const PlaneInfo &    getPlaneInfo(const unsigned int d) const;
    const MPI_Comm &     getCommGrid(void) const;
    const MPI_Comm &     getDirCommGrid(const unsigned int d) const;
    RngType &            getRng(const unsigned int i) const;
    // check
    static void checkDir(const unsigned int d);
    // get absolute direction
    static unsigned int absDir(const unsigned int d);
    // get opposite direction
    static unsigned int oppDir(const unsigned int d);
    // get neighbor coordinate in the grid
    int neighborCoord(const unsigned int d) const;
    // need communication in direction d?
    bool needComm(const unsigned d) const;
    // clock
    static double time(void);
    // random generator management
    void generateSeed(void);
    void initializeRng(void);
private:
    int                                              rank_;
    unsigned int                                     volume_, locVolume_;
    unsigned int                                     commBufferSize_;
    std::array<int, D>                               p_, coord_;
    std::array<bool, D>                              needComm_;
    Coord<D>                                         locSurface_, dim_, locDim_;
    MPI_Comm                                         commGrid_;
    std::array<MPI_Comm, D>                          dirCommGrid_;
    std::array<PlaneInfo, D>                         planeInfo_;
    std::array<std::unique_ptr<unsigned int[]>, 2*D> nnTable_;
    std::unique_ptr<SeedType[]>                      seed_;
    std::unique_ptr<RngType[]>                       rng_;
};

// global layout
extern LayoutObject *globalLayout;

template <unsigned int D>
void registerGlobalLayout(Layout<D> &layout)
{
    globalLayout = &layout;
}

/******************************************************************************
 *                           Layout implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <unsigned int D>
Layout<D>::Layout(const Coord<D> &dim, const Coord<D> &p)
: Logger("Layout")
, dim_(dim)
{
    int          size, isInit, nproc = 1, isPeriodic[D], isActive[D];
    unsigned int shift = 0;
    std::string  buf;

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
        if (p[d] == 0)
        {
            locGlobalError("no processes in direction " + strFrom(d));
        }
        p_[d]        = static_cast<int>(p[d]);
        needComm_[d] = (p_[d] != 1);
    }
    buf = "";
    for (unsigned int d = 0; d < D; ++d)
    {
        buf += strFrom(dim_[d]) + ((d != D - 1) ? "x" : "");
    }
    masterLog("lattice size     : " + buf);
    buf = "";
    for (unsigned int d = 0; d < D; ++d)
    {
        buf += strFrom(p_[d]) + ((d != D - 1) ? "x" : "");
    }
    masterLog("MPI partition    : " + buf);
    for (unsigned int d = 0; d < D; ++d)
    {
        if (dim_[d] % static_cast<unsigned int>(p_[d]) != 0)
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
    volume_         = 1;
    locVolume_      = 1;
    commBufferSize_ = 0;
    for (unsigned int d = 0; d < D; ++d)
    {
        locDim_[d] = dim_[d]/static_cast<unsigned int>(p_[d]);
    }
    for (unsigned int d = 0; d < D; ++d)
    {
        volume_        *= dim_[d];
        locVolume_     *= locDim_[d];
        locSurface_[d]  = 1;
        for (unsigned int i = 0; i < D; ++i)
        {
            locSurface_[d] *= (i != d) ? locDim_[i] : 1;
        }
        commBufferSize_ += 2*locSurface_[d];
    }
    masterLog("volume           : " + strFrom(volume_));
    masterLog("local volume     : " + strFrom(locVolume_));
    buf = "";
    for (unsigned int d = 0; d < D; ++d)
    {
        buf += (d == 0) ? "[" : "";
        buf += strFrom(locSurface_[d]);
        buf += (d == D - 1) ? "]" : ", ";
    }
    masterLog("local surfaces   : " + buf);
    masterLog("comm buffer size : " + strFrom(commBufferSize_));
    masterLog("total local size : " + strFrom(locVolume_ + commBufferSize_));

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
        MPI_Cart_sub(commGrid_, isActive, &dirCommGrid_[d]);
    }
    MPI_Comm_rank(commGrid_, &rank_);
    MPI_Cart_coords(commGrid_, rank_, D, coord_.data());

    // compute plane communication informations
    planeInfo_[0].nBlocks   = 1;
    planeInfo_[0].blockSize = 1;
    planeInfo_[0].stride    = 1;
    for (unsigned int d = 1; d < D; ++d)
    {
        planeInfo_[0].blockSize *= locDim_[d];
    }
    for (unsigned int d = 1; d < D; ++d)
    {
        planeInfo_[d].nBlocks   = planeInfo_[d-1].nBlocks*locDim_[d-1];
        planeInfo_[d].blockSize = planeInfo_[d-1].blockSize/locDim_[d];
        planeInfo_[d].stride    = planeInfo_[d-1].blockSize;
    }
    masterLog("plane indexing   :");
    for (unsigned int d = 0; d < D; ++d)
    {
        masterLog("  dir= " + strFrom(d) + ": nBlocks= " +
                  strFrom(planeInfo_[d].nBlocks) + ", size= " +
                  strFrom(planeInfo_[d].blockSize) + ", stride= " +
                  strFrom(planeInfo_[d].stride));
    }

    // make nearest-neighbour table
    shift = locVolume_;
    for (unsigned int d = 0; d < 2*D; ++d)
    {
        Coord<D>           x;
        Coord<D-1>         planeX, planeDim;
        const unsigned int ad = absDir(d);

        nnTable_[d].reset(new unsigned int[locVolume_]);
        for (unsigned int e = 0; e < D - 1; ++e)
        {
            planeDim[e] = (e < ad) ? locDim_[e] : locDim_[e+1];
        }
        for (unsigned int i = 0; i < locVolume_; ++i)
        {
            x = rowMajorToCoord(i, locDim_);
            // positive direction, local or no communications
            if ((ad == d)&&((x[ad] != locDim_[ad] - 1)||!needComm_[ad]))
            {
                x[ad]          = (x[ad] + 1 + locDim_[ad]) % locDim_[ad];
                nnTable_[d][i] = coordToRowMajor(x, locDim_);
            }
            // negative direction, local or no communications
            else if ((ad != d)&&((x[ad] != 0)||!needComm_[ad]))
            {
                x[ad]          = (x[ad] - 1 + locDim_[ad]) % locDim_[ad];
                nnTable_[d][i] = coordToRowMajor(x, locDim_);
            }
            // mapping to communication buffer
            else
            {
                for (unsigned int e = 0; e < D - 1; ++e)
                {
                    planeX[e] = (e < ad) ? x[e] : x[e+1];
                }
                nnTable_[d][i]  = coordToRowMajor(planeX, planeDim);
                nnTable_[d][i] += shift;
            }

        }
        shift += locSurface_[ad];
    }

    // make random generator lattice
    seed_.reset(new SeedType[locVolume_]);
    rng_.reset(new RngType[locVolume_]);
    generateSeed();
    initializeRng();

    // everything is ok
    masterLog("communication grid ready");
}

// destructor //////////////////////////////////////////////////////////////////
template <unsigned int D>
Layout<D>::~Layout(void)
{
    MPI_Barrier(MPI_COMM_WORLD);
    masterLog("communication grid terminated");
    MPI_Finalize();
}

// access //////////////////////////////////////////////////////////////////////
template <unsigned int D>
unsigned int Layout<D>::getRank(void) const
{
    return static_cast<unsigned int>(rank_);
}

template <unsigned int D>
const Coord<D> & Layout<D>::getDim(void) const
{
    return dim_;
}

template <unsigned int D>
unsigned int Layout<D>::getDim(const unsigned int d) const
{
    checkDir(d);

    return dim_[absDir(d)];
}

template <unsigned int D>
const Coord<D> & Layout<D>::getLocalDim(void) const
{
    return locDim_;
}

template <unsigned int D>
unsigned int Layout<D>::getLocalDim(const unsigned int d) const
{
    checkDir(d);

    return locDim_[absDir(d)];
}

template <unsigned int D>
const Coord<D> & Layout<D>::getLocalSurface(void) const
{
    return locSurface_;
}

template <unsigned int D>
unsigned int Layout<D>::getLocalSurface(const unsigned int d) const
{
    checkDir(d);

    return locSurface_[absDir(d)];
}

template <unsigned int D>
unsigned int Layout<D>::getVolume(void) const
{
    return volume_;
}

template <unsigned int D>
unsigned int Layout<D>::getLocalVolume(void) const
{
    return locVolume_;
}

template <unsigned int D>
unsigned int Layout<D>::getCommBufferSize(void) const
{
    return commBufferSize_;
}

template <unsigned int D>
unsigned int Layout<D>::getNearNeigh(const unsigned int i,
                                     const unsigned int d) const
{
    return nnTable_[d][i];
}

template <unsigned int D>
const unsigned int * Layout<D>::getNearNeigh(const unsigned int d) const
{
    return nnTable_[d].get();
}

template <unsigned int D>
const typename Layout<D>::PlaneInfo &
Layout<D>::getPlaneInfo(const unsigned int d) const
{
    checkDir(d);

    return planeInfo_[absDir(d)];
}

template <unsigned int D>
const MPI_Comm & Layout<D>::getCommGrid(void) const
{
    return commGrid_;
}

template <unsigned int D>
const MPI_Comm & Layout<D>::getDirCommGrid(const unsigned int d) const
{
    checkDir(d);

    return dirCommGrid_[absDir(d)];
}

template <unsigned int D>
typename Layout<D>::RngType & Layout<D>::getRng(const unsigned int i) const
{
    return rng_[i];
}

// check ///////////////////////////////////////////////////////////////////////
template <unsigned int D>
void Layout<D>::checkDir(const unsigned int d)
{
    if (d >= 2*D)
    {
        locGlobalError("direction " + strFrom(d) + " invalid " +
                       "(number of dimension is " + strFrom(D) + ")");
    }
}

// get absolute direction //////////////////////////////////////////////////////
template <unsigned int D>
unsigned int Layout<D>::absDir(const unsigned int d)
{
    return d % D;
}

// get opposite direction //////////////////////////////////////////////////////
template <unsigned int D>
unsigned int Layout<D>::oppDir(const unsigned int d)
{
    return (d + D) % (2*D);
}

// get neighbor coordinate in the grid /////////////////////////////////////////
template <unsigned int D>
int Layout<D>::neighborCoord(const unsigned int d) const
{
    unsigned int ad = absDir(d);

    checkDir(d);
    
    return (d < D) ? ((coord_[ad] + 1)          % p_[ad])
                   : ((coord_[ad] + p_[ad] - 1) % p_[ad]);
}

// need communication in direction d? //////////////////////////////////////////
template <unsigned int D>
bool Layout<D>::needComm(const unsigned d) const
{
    return needComm_[absDir(d)];
}

// clock ///////////////////////////////////////////////////////////////////////
template <unsigned int D>
double Layout<D>::time(void)
{
    return MPI_Wtime();
}

// random generator management /////////////////////////////////////////////////
template <unsigned int D>
void Layout<D>::generateSeed(void)
{
    std::random_device rDev;

    for (unsigned int i = 0; i < locVolume_; ++i)
    {
        seed_[i] = static_cast<SeedType>(rDev());
    }
}

template <unsigned int D>
void Layout<D>::initializeRng(void)
{
    for (unsigned int i = 0; i < locVolume_; ++i)
    {
        rng_[i].seed(seed_[i]);
    }
}

END_NAMESPACE

#endif // LatSim_Layout_hpp_
