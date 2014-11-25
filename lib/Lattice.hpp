/*
 * Lattice.hpp, part of LatSim
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

#ifndef LatSim_Lattice_hpp_
#define LatSim_Lattice_hpp_

#include <LatSim/Global.hpp>
#include <LatSim/Expressions.hpp>
#include <LatSim/Layout.hpp>
#include <LatSim/Logger.hpp>
#include <LatSim/MpiTypes.hpp>

BEGIN_NAMESPACE

/******************************************************************************
 *                                Lattice                                     *
 ******************************************************************************/

template <typename T, unsigned long D>
class Lattice: public Logger
{
public:
    typedef T SiteType;
public:
    // constructors
    explicit Lattice(const LayoutObject *layout = globalLayout);
    Lattice(const Lattice<T, D> &l);
    // destructor
    virtual ~Lattice(void);
    // local site access
    inline const T & operator[](const unsigned long i) const;
    inline T &       operator[](const unsigned long i);
    // directional gathering
    //// non-blocking
    void gatherStart(const unsigned int d);
    void gatherWait(const unsigned int d);
    //// blocking
    void gather(const unsigned d);
    // assignement operator
    template <typename... Ts, typename... Ops>
    Lattice<T, D> & operator=(const std::pair<std::tuple<Ts...>,
                                              std::tuple<Ops...>>& stack);
private:
    // helpers for constructors
    void reallocate(const LayoutObject *layout = globalLayout);
    void createMpiTypes(void);
private:
    const Layout<D>              *layout_;
    std::unique_ptr<T>           data_;
    T                            *lattice_;
    std::array<T *, 2*D>         commBuffer_;
    MPI_Datatype                 mpiElemType_;
    std::array<MPI_Datatype, D>  mpiBufType_, mpiPlaneType_;
    std::array<MPI_Request, 2*D> sReq_, rReq_;
    std::array<MPI_Status, 2*D>  sStatus_, rStatus_;
};

// arithmetic operators
#define DEFINE_OP(op, lname, rname)\
template <typename T, unsigned long D>\
std::pair<std::tuple<const Lattice<T, D> &, const Lattice<T, D> &>,\
          std::tuple<Expr::rname<T, T>>>\
op(const Lattice<T, D> &lhs, const Lattice<T, D> &rhs)\
{\
    return std::make_pair(\
        std::tuple<const Lattice<T, D> &, const Lattice<T, D> &>(lhs, rhs),\
        std::tuple<Expr::rname<T, T>>(Expr::rname<T, T>()));\
}\
\
template <typename... Ts, typename... Ops, typename T, unsigned long D>\
std::pair<std::tuple<Ts..., const Lattice<T, D> &>,\
          std::tuple<Ops..., Expr::rname<T, T>>>\
op(const std::pair<std::tuple<Ts...>, std::tuple<Ops...>> &lhs,\
   const Lattice<T, D> &rhs)\
{\
    return std::make_pair(\
        std::tuple_cat(lhs.first, std::tuple<const Lattice<T, D> &>(rhs)),\
        std::tuple_cat(lhs.second,\
                       std::tuple<Expr::rname<T, T>>(Expr::rname<T, T>())));\
}\
\
template <typename... Ts, typename... Ops, typename T, unsigned long D>\
std::pair<std::tuple<Ts..., const Lattice<T, D> &>,\
          std::tuple<Ops..., Expr::lname<T, T>>>\
op(const Lattice<T, D> &lhs,\
   const std::pair<std::tuple<Ts...>, std::tuple<Ops...>> &rhs)\
{\
    return std::make_pair(\
        std::tuple_cat(rhs.first, std::tuple<const Lattice<T, D> &>(lhs)),\
        std::tuple_cat(rhs.second,\
                       std::tuple<Expr::lname<T, T>>(Expr::lname<T, T>())));\
}

DEFINE_OP(operator+, Add,  Add)
DEFINE_OP(operator-, LSub, RSub)
DEFINE_OP(operator*, LMul, RMul)
DEFINE_OP(operator/, LDiv, RDiv)

/******************************************************************************
 *                          Lattice implementation                            *
 ******************************************************************************/
// helpers for constructors ////////////////////////////////////////////////////
template <typename T, unsigned long D>
void Lattice<T, D>::reallocate(const LayoutObject *layout)
{
    // allocate lattice and communication buffers
    layout_ = dynamic_cast<const Layout<D> *>(layout);
    data_.reset(new T[layout_->getLocalVolume()+layout_->getCommBufferSize()]);
    lattice_       = data_.get();
    commBuffer_[0] = data_.get() + layout_->getLocalVolume();
    for (unsigned int d = 1; d < 2*D; ++d)
    {
        commBuffer_[d] = commBuffer_[d-1] + layout_->getLocalSurface(d);
    }
}

template <typename T, unsigned long D>
void Lattice<T, D>::createMpiTypes(void)
{
    // create base MPI data type
    MpiType<T>::make(mpiElemType_);
    MPI_Type_commit(&mpiElemType_);

    // creat MPI types for planes and communication buffers
    for (unsigned int d = 0; d < D; ++d)
    {
        const auto &p = layout_->getPlaneInfo(d);

        MPI_Type_contiguous(layout_->getLocalSurface(d), mpiElemType_,
                            &mpiBufType_[d]);
        MPI_Type_commit(&mpiBufType_[d]);
        MPI_Type_vector(p.nBlocks, p.blockSize, p.stride, mpiElemType_,
                        &mpiPlaneType_[d]);
        MPI_Type_commit(&mpiPlaneType_[d]);
    }
}

// constructors ////////////////////////////////////////////////////////////////
template <typename T, unsigned long D>
Lattice<T, D>::Lattice(const LayoutObject *layout)
: Logger("Lattice")
{
    reallocate(layout);
    createMpiTypes();
}

template <typename T, unsigned long D>
Lattice<T, D>::Lattice(const Lattice<T, D> &rhs)
{
    if (this != &rhs)
    {
        const unsigned long size = rhs.layout_->getLocalVolume();

        if (layout_ != rhs.layout_)
        {
            reallocate(rhs.layout_);
        }
        std::copy(rhs.lattice_, rhs.lattice_ + size, size);
        createMpiTypes();
    }
}

// destructor //////////////////////////////////////////////////////////////////
template <typename T, unsigned long D>
Lattice<T, D>::~Lattice(void)
{
    for (unsigned int d = 0; d < D; ++d)
    {
        MPI_Type_free(&mpiBufType_[d]);
        MPI_Type_free(&mpiPlaneType_[d]);
    }
    MPI_Type_free(&mpiElemType_);
}

// local site access
template <typename T, unsigned long D>
const T & Lattice<T, D>::operator[](const unsigned long i) const
{
    if (i >= layout_->getLocalVolume())
    {
        locGlobalError("site index out of range (i= " + strFrom(i) + ")");
    }
    
    return lattice_[i];
}

template <typename T, unsigned long D>
T & Lattice<T, D>::operator[](const unsigned long i)
{
    if (i >= layout_->getLocalVolume())
    {
        locGlobalError("site index out of range (i= " + strFrom(i) + ")");
    }
    
    return lattice_[i];
}

// directional gathering ///////////////////////////////////////////////////////
template <typename T, unsigned long D>
void Lattice<T, D>::gatherStart(const unsigned int d)
{
    if (layout_->needComm(d))
    {
        const unsigned int pShift  = (layout_->getLocalDim(d) - 1)
                                     *layout_->getPlaneInfo(d).blockSize;
        const int          ad      = layout_->absDir(d);
        const T            *sendPt = lattice_ + ((d < D) ? 0 : pShift);

        MPI_Isend(sendPt, 1, mpiPlaneType_[ad],
                  layout_->neighborCoord(layout_->oppDir(d)), d,
                  layout_->getDirCommGrid(d), &sReq_[d]);
        MPI_Irecv(commBuffer_[d], 1, mpiBufType_[ad],
                  layout_->neighborCoord(d), d, layout_->getDirCommGrid(d),
                  &rReq_[d]);
    }
}

template <typename T, unsigned long D>
void Lattice<T, D>::gatherWait(const unsigned int d)
{
    if (layout_->needComm(d))
    {
        MPI_Wait(&sReq_[d], &sStatus_[d]);
        MPI_Wait(&rReq_[d], &rStatus_[d]);
    }
}

template <typename T, unsigned long D>
void Lattice<T, D>::gather(const unsigned int d)
{
    gatherStart(d);
    gatherWait(d);
}

// assignement operator ////////////////////////////////////////////////////////
template <typename T, unsigned long D>
template <typename... Ts, typename... Ops>
Lattice<T, D> &
Lattice<T, D>::operator=(const std::pair<std::tuple<Ts...>,
                                         std::tuple<Ops...>>& stack)
{
    for (unsigned int i = 0; i < layout_->getLocalVolume(); ++i)
    {
        lattice_[i] = Expr::eval<T>(Expr::getSites(i, stack.first),
                                    stack.second);
    }
    
    return *this;
}

END_NAMESPACE

#endif // LatSim_Lattice_hpp_
