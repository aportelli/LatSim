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
    template <typename Op, typename... Ts>
    Lattice<T, D> & operator=(const std::tuple<Op, Ts...> &expr) flatten
    {
        for (unsigned long i = 0; i < layout_->getLocalVolume(); ++i)
        {
            lattice_[i] = Expr::eval(i, expr);
        }

        return *this;
    }
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
#define OP_NAME(name) Expr::name<decltype(Expr::eval(0lu, lhs)),\
                                 decltype(Expr::eval(0lu, rhs))>

#define DEFINE_OP(op, name)\
template <typename T1, typename T2, unsigned long D>\
strong_inline auto op(const Lattice<T1, D> &lhs, const Lattice<T2, D> &rhs)\
->std::tuple<OP_NAME(name), const Lattice<T1, D> &,\
             const Lattice<T2, D> &>\
{\
    return std::tuple<OP_NAME(name), const Lattice<T1, D> &,\
        const Lattice<T2, D> &>(OP_NAME(name)(), lhs, rhs);\
}\
template <typename Op1, typename... T1s, typename T2, unsigned long D>\
strong_inline auto op(const std::tuple<Op1, T1s...> &lhs,\
                     const Lattice<T2, D> &rhs)\
->std::tuple<OP_NAME(name), const std::tuple<Op1, T1s...> &,\
             const Lattice<T2, D> &>\
{\
    return std::tuple<OP_NAME(name), const std::tuple<Op1, T1s...> &,\
        const Lattice<T2, D> &>(OP_NAME(name)(), lhs, rhs);\
}\
template <typename T1, typename Op2, typename... T2s, unsigned long D>\
strong_inline auto op(const Lattice<T1, D> &lhs,\
                     const std::tuple<Op2, T2s...> &rhs)\
->std::tuple<OP_NAME(name), const Lattice<T1, D> &,\
             const std::tuple<Op2, T2s...> &>\
{\
    return std::tuple<OP_NAME(name), const Lattice<T1, D> &,\
        const std::tuple<Op2, T2s...> &>(OP_NAME(name)(), lhs, rhs);\
}\
template <typename Op1, typename... T1s, typename Op2, typename... T2s>\
strong_inline auto op(const std::tuple<Op1, T1s...> &lhs,\
               const std::tuple<Op2, T2s...> &rhs)\
->std::tuple<OP_NAME(name), const std::tuple<Op1, T1s...> &,\
             const std::tuple<Op2, T2s...> &>\
{\
    return std::tuple<OP_NAME(name), const std::tuple<Op1, T1s...> &,\
        const std::tuple<Op2, T2s...> &>(OP_NAME(name)(), lhs, rhs);\
}

DEFINE_OP(operator+, Add);
DEFINE_OP(operator-, Sub);
DEFINE_OP(operator*, Mul);
DEFINE_OP(operator/, Div);

/******************************************************************************
 *                          Lattice implementation                            *
 ******************************************************************************/
// helpers for constructors ////////////////////////////////////////////////////
template <typename T, unsigned long D>
void Lattice<T, D>::reallocate(const LayoutObject *layout)
{
    // set layout
    layout_  = dynamic_cast<const Layout<D> *>(layout);

    // allocate lattice and communication buffers
    unsigned long size = layout_->getLocalVolume()+layout_->getCommBufferSize();
    void          *start;

    size += alignof(T);
    data_.reset(new T[size]);
    start = data_.get();
    std::align(alignof(T), size*sizeof(T), start, size);
    lattice_ = static_cast<T *>(start);
    commBuffer_[0] = lattice_ + layout_->getLocalVolume();
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
inline const T & Lattice<T, D>::operator[](const unsigned long i) const
{
    return lattice_[i];
}

template <typename T, unsigned long D>
inline T & Lattice<T, D>::operator[](const unsigned long i)
{
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

END_NAMESPACE

#endif // LatSim_Lattice_hpp_
