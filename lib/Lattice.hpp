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
template <typename T, unsigned int D>
class Lattice: public Logger, public LatticeObj
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
    inline const T & operator[](const unsigned int i) const;
    inline T &       operator[](const unsigned int i);
    // directional gathering
    //// non-blocking
    void gatherStart(const unsigned int d);
    void gatherWait(const unsigned int d);
    //// blocking
    void gather(const unsigned int d);
    // expression evaluation
    template <typename Op, typename... Ts>
    inline Lattice<T, D> & operator=(const LatExpr<Op, Ts...> &expr) flatten;
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
#define OP(name) Expr::name<decltype(Expr::eval(0lu, lhs)),\
                                 decltype(Expr::eval(0lu, rhs))>
#define DEFINE_OP(op, name)\
template <typename T1, typename T2>\
strong_inline auto op(const T1 &lhs, const T2 &rhs)\
->typename std::enable_if<std::is_base_of<LatticeObj, T1>::value\
                          &&std::is_base_of<LatticeObj, T2>::value,\
                          LatExpr<OP(name), const T1 &, const T2 &>>::type\
{\
    return LatExpr<OP(name), const T1 &, const T2 &>(OP(name)(), lhs, rhs);\
}

DEFINE_OP(operator+, Add)
DEFINE_OP(operator-, Sub)
DEFINE_OP(operator*, Mul)
DEFINE_OP(operator/, Div)

/******************************************************************************
 *                          Lattice implementation                            *
 ******************************************************************************/
// helpers for constructors ////////////////////////////////////////////////////
template <typename T, unsigned int D>
void Lattice<T, D>::reallocate(const LayoutObject *layout)
{
    // set layout
    layout_  = dynamic_cast<const Layout<D> *>(layout);

    // allocate lattice and communication buffers
    unsigned int size = layout_->getLocalVolume()+layout_->getCommBufferSize();

    size += alignof(T);
    data_.reset(new T[size]);
    lattice_       = data_.get();
    commBuffer_[0] = lattice_ + layout_->getLocalVolume();
    for (unsigned int d = 1; d < 2*D; ++d)
    {
        commBuffer_[d] = commBuffer_[d-1] + layout_->getLocalSurface(d);
    }
}

template <typename T, unsigned int D>
void Lattice<T, D>::createMpiTypes(void)
{
    // create base MPI data type
    MpiType<T>::make(mpiElemType_);
    MPI_Type_commit(&mpiElemType_);

    // creat MPI types for planes and communication buffers
    for (unsigned int d = 0; d < D; ++d)
    {
        const auto &p         = layout_->getPlaneInfo(d);
        const int  locSurface = static_cast<int>(layout_->getLocalSurface(d));
        const int  n          = static_cast<int>(p.nBlocks);
        const int  size       = static_cast<int>(p.blockSize);
        const int  stride     = static_cast<int>(p.stride);

        MPI_Type_contiguous(locSurface, mpiElemType_, &mpiBufType_[d]);
        MPI_Type_commit(&mpiBufType_[d]);
        MPI_Type_vector(n, size, stride, mpiElemType_, &mpiPlaneType_[d]);
        MPI_Type_commit(&mpiPlaneType_[d]);
    }
}

// constructors ////////////////////////////////////////////////////////////////
template <typename T, unsigned int D>
Lattice<T, D>::Lattice(const LayoutObject *layout)
: Logger("Lattice")
{
    reallocate(layout);
    createMpiTypes();
}

template <typename T, unsigned int D>
Lattice<T, D>::Lattice(const Lattice<T, D> &rhs)
{
    if (this != &rhs)
    {
        const unsigned int size = rhs.layout_->getLocalVolume();

        if (layout_ != rhs.layout_)
        {
            reallocate(rhs.layout_);
        }
        std::copy(rhs.lattice_, rhs.lattice_ + size, size);
        createMpiTypes();
    }
}

// destructor //////////////////////////////////////////////////////////////////
template <typename T, unsigned int D>
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
template <typename T, unsigned int D>
inline const T & Lattice<T, D>::operator[](const unsigned int i) const
{
    return lattice_[i];
}

template <typename T, unsigned int D>
inline T & Lattice<T, D>::operator[](const unsigned int i)
{
    return lattice_[i];
}

// directional gathering ///////////////////////////////////////////////////////
template <typename T, unsigned int D>
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

template <typename T, unsigned int D>
void Lattice<T, D>::gatherWait(const unsigned int d)
{
    if (layout_->needComm(d))
    {
        MPI_Wait(&sReq_[d], &sStatus_[d]);
        MPI_Wait(&rReq_[d], &rStatus_[d]);
    }
}

template <typename T, unsigned int D>
void Lattice<T, D>::gather(const unsigned int d)
{
    gatherStart(d);
    gatherWait(d);
}

// expression evaluation ///////////////////////////////////////////////////////
template <typename T, unsigned int D>
template <typename Op, typename... Ts>
inline Lattice<T, D> &
Lattice<T, D>::operator=(const LatExpr<Op, Ts...> &expr) flatten
{
    for (unsigned int i = 0; i < layout_->getLocalVolume(); ++i)
    {
        lattice_[i] = Expr::eval(i, expr);
    }

    return *this;
}

END_NAMESPACE

#endif // LatSim_Lattice_hpp_
