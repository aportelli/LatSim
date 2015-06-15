/*
 * Lattice.hpp, part of LatSim
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

#ifndef LatSim_Lattice_hpp_
#define LatSim_Lattice_hpp_

#include <LatSim/Global.hpp>
#include <LatSim/Expressions.hpp>
#include <LatSim/Layout.hpp>
#include <LatSim/Logger.hpp>
#include <LatSim/MpiTypes.hpp>

BEGIN_LATSIM_NAMESPACE

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
    Lattice(Lattice<T, D> &&l);
    // destructor
    virtual ~Lattice(void);
    // local site access
    inline const T & operator()(const unsigned int i) const;
    inline T &       operator()(const unsigned int i);
    inline const T & operator()(const unsigned int i,
                                const unsigned int d) const;
    inline T &       operator()(const unsigned int i, const unsigned int d);
    // global site access
    T getSite(const Coord<D> &x);
    // filling function
    void fill(const T &x);
    // layout access
    const Layout<D> & getLayout(void);
    // directional gathering
    //// non-blocking
    void gatherStart(const unsigned int d);
    void gatherWait(const unsigned int d);
    //// blocking
    void gather(const unsigned int d);
    // assignement operator
    Lattice<T, D> & operator=(const Lattice<T, D> &l);
    // expression evaluation
    template <typename Op, typename... Ts>
    strong_inline flatten Lattice<T, D> & operator=(const LatExpr<Op, Ts...> &expr);
    template <typename Op>
    strong_inline flatten Lattice<T, D> &
    operator=(const LatOpExpr<Op, Lattice<T, D>> &expr) ;
    // dump to stdout
    void dump(void);
    // reduction
    T reduce(const T &x, std::function<T(const Vec<T> &buf)> &f);
    T reduceSum(const T &x);
    // lattice integration
    T      sum(void);
    Vec<T> sumSlice(const unsigned int d);
private:
    // helpers for constructors/destructor
    void reallocate(const LayoutObject *layout = globalLayout);
    void copy(const Lattice<T, D> &l);
    void clear(void);
    void swap(Lattice<T, D> &l);
    void createMpiTypes(void);
private:
    const Layout<D>              *layout_;
    std::unique_ptr<T[]>         data_;
    T                            *lattice_;
    std::array<T *, 2*D>         commBuffer_;
    Vec<T>                       reduceBuffer_;
    bool                         mpiTypesInit_;
    MPI_Datatype                 mpiSiteType_;
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

// useful loop macro
#define FOR_SITE(l, i)\
for (unsigned int i = 0; i < (l).getLayout().getLocalVolume(); ++i)

#define PINFO(l,d) (l).getLayout().getPlaneInfo(d)
#define FOR_SITE_IN_PLANE(l, d, x, i)\
for (unsigned int _b = 0; _b < PINFO(l,d).nBlocks; ++_b)\
for (unsigned int i = (x)*PINFO(l,d).blockSize + PINFO(l,d).stride*_b;\
     i < (x+1)*PINFO(l,d).blockSize + PINFO(l,d).stride*_b; ++i)

/******************************************************************************
 *                          Lattice implementation                            *
 ******************************************************************************/
// helpers for constructors ////////////////////////////////////////////////////
template <typename T, unsigned int D>
void Lattice<T, D>::reallocate(const LayoutObject *layout)
{
    // set layout
    layout_ = dynamic_cast<const Layout<D> *>(layout);

    // allocate lattice and communication buffers
    unsigned int size = layout_->getLocalVolume()+layout_->getCommBufferSize();

    size += alignof(T);
    data_.reset(new T[size]);
    lattice_       = data_.get();
    commBuffer_[0] = lattice_ + layout_->getLocalVolume();
    for (unsigned int d = 1; d < 2*D; ++d)
    {
        commBuffer_[d] = commBuffer_[d-1] + layout_->getLocalSurface(d-1);
    }
    reduceBuffer_.resize(static_cast<Index>(layout_->getNProcess()));
}

template <typename T, unsigned int D>
void Lattice<T, D>::clear(void)
{
    layout_  = nullptr;
    data_.reset(nullptr);
    lattice_ = nullptr;
    for (unsigned int d = 0; d < D; ++d)
    {
        commBuffer_[d] = nullptr;
    }
    if (mpiTypesInit_)
    {
        for (unsigned int d = 0; d < D; ++d)
        {
            MPI_Type_free(&mpiBufType_[d]);
            MPI_Type_free(&mpiPlaneType_[d]);
        }
        MPI_Type_free(&mpiSiteType_);
    }
    mpiTypesInit_ = false;
}

template <typename T, unsigned int D>
void Lattice<T, D>::copy(const Lattice<T, D> &rhs)
{
    if (this != &rhs)
    {
        const unsigned int size = rhs.layout_->getLocalVolume();

        if (layout_ != rhs.layout_)
        {
            reallocate(rhs.layout_);
            createMpiTypes();
        }
        std::copy(rhs.lattice_, rhs.lattice_ + size, lattice_);
    }
}

template <typename T, unsigned int D>
void Lattice<T, D>::swap(Lattice<T, D> &l)
{
    using std::swap;

    swap(layout_, l.layout_);
    swap(data_, l.data_);
    swap(lattice_, l.lattice_);
    swap(commBuffer_, l.commBuffer_);
    reduceBuffer_.swap(l.reduceBuffer_);
    swap(mpiTypesInit_, l.mpiTypesInit_);
    swap(mpiSiteType_, l.mpiSiteType_);
    swap(mpiBufType_, l.mpiBufType_);
    swap(mpiPlaneType_, l.mpiPlaneType_);
    swap(sReq_, l.sReq_);
    swap(rReq_, l.rReq_);
    swap(sStatus_, l.sStatus_);
    swap(rStatus_, l.rStatus_);
}

template <typename T, unsigned int D>
void Lattice<T, D>::createMpiTypes(void)
{
    // create base MPI data type
    MpiType<T>::make(mpiSiteType_);
    MPI_Type_commit(&mpiSiteType_);

    // create MPI types for planes and communication buffers
    for (unsigned int d = 0; d < D; ++d)
    {
        const auto &p         = layout_->getPlaneInfo(d);
        const int  locSurface = static_cast<int>(layout_->getLocalSurface(d));
        const int  n          = static_cast<int>(p.nBlocks);
        const int  size       = static_cast<int>(p.blockSize);
        const int  stride     = static_cast<int>(p.stride);

        MPI_Type_contiguous(locSurface, mpiSiteType_, &mpiBufType_[d]);
        MPI_Type_commit(&mpiBufType_[d]);
        MPI_Type_vector(n, size, stride, mpiSiteType_, &mpiPlaneType_[d]);
        MPI_Type_commit(&mpiPlaneType_[d]);
    }
    mpiTypesInit_ = true;
}

// constructors ////////////////////////////////////////////////////////////////
template <typename T, unsigned int D>
Lattice<T, D>::Lattice(const LayoutObject *layout)
: Logger("Lattice")
{
    if (layout != nullptr)
    {
        reallocate(layout);
        createMpiTypes();
    }
    else
    {
        clear();
    }
}

template <typename T, unsigned int D>
Lattice<T, D>::Lattice(const Lattice<T, D> &rhs)
: Lattice(nullptr)
{
    copy(rhs);
}

template <typename T, unsigned int D>
Lattice<T, D>::Lattice(Lattice<T, D> &&rhs)
: Lattice(nullptr)
{
    this->swap(rhs);
}

// destructor //////////////////////////////////////////////////////////////////
template <typename T, unsigned int D>
Lattice<T, D>::~Lattice(void)
{
    clear();
}

// global site access //////////////////////////////////////////////////////////
template <typename T, unsigned int D>
T Lattice<T, D>::getSite(const Coord<D> &x)
{
    unsigned int siteRank = layout_->getSiteRank(x);
    T            site;

    if (layout_->getMyRank() == siteRank)
    {
        Coord<D> locX;

        for (unsigned int d = 0; d < D; ++d)
        {
            locX[d] = x[d] - layout_->getFirstSite()[d];
        }
        site = (*this)(layout_->getIndex(locX));
    }
    MPI_Bcast(&site, 1, mpiSiteType_, siteRank, layout_->getCommGrid());

    return site;
}

// local site access ///////////////////////////////////////////////////////////
template <typename T, unsigned int D>
inline const T & Lattice<T, D>::operator()(const unsigned int i) const
{
    return lattice_[i];
}

template <typename T, unsigned int D>
inline T & Lattice<T, D>::operator()(const unsigned int i)
{
    return lattice_[i];
}

template <typename T, unsigned int D>
inline const T & Lattice<T, D>::operator()(const unsigned int i,
                                           const unsigned int d) const
{
    return lattice_[layout_->getNearNeigh(i, d)];
}

template <typename T, unsigned int D>
inline T & Lattice<T, D>::operator()(const unsigned int i, const unsigned int d)
{
    return lattice_[layout_->getNearNeigh(i, d)];
}

// filling function ////////////////////////////////////////////////////////////
template <typename T, unsigned int D>
void Lattice<T, D>::fill(const T &x)
{
    FOR_SITE(*this, i)
    {
        (*this)(i) = x;
    }
}

// layout access ///////////////////////////////////////////////////////////////
template <typename T, unsigned int D>
const Layout<D> & Lattice<T, D>::getLayout(void)
{
    return *layout_;
}

// directional gathering ///////////////////////////////////////////////////////
template <typename T, unsigned int D>
void Lattice<T, D>::gatherStart(const unsigned int d)
{
    const unsigned int pShift  = (layout_->getLocalDim(d) - 1)
                                 *layout_->getPlaneInfo(d).blockSize;
    const int sd               = static_cast<int>(d);
    const unsigned int ad      = layout_->absDir(d);
    T                  *sendPt = lattice_ + ((d < D) ? 0 : pShift);

    MPI_Isend(sendPt, 1, mpiPlaneType_[ad],
              layout_->neighborCoord(layout_->oppDir(d)), sd,
              layout_->getDirCommGrid(d), &sReq_[d]);
    MPI_Irecv(commBuffer_[d], 1, mpiBufType_[ad],
              layout_->neighborCoord(d), sd, layout_->getDirCommGrid(d),
              &rReq_[d]);
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

// assignement operator ////////////////////////////////////////////////////////
template <typename T, unsigned int D>
Lattice<T, D> & Lattice<T, D>::operator=(const Lattice<T, D> &l)
{
    copy(l);

    return *this;
}

// expression evaluation ///////////////////////////////////////////////////////
template <typename T, unsigned int D>
template <typename Op, typename... Ts>
strong_inline flatten Lattice<T, D> &
Lattice<T, D>::operator=(const LatExpr<Op, Ts...> &expr)
{
    for (unsigned int i = 0; i < layout_->getLocalVolume(); ++i)
    {
        lattice_[i] = Expr::eval(i, expr);
    }

    return *this;
}

template <typename T, unsigned int D>
template <typename Op>
strong_inline flatten Lattice<T, D> &
Lattice<T, D>::operator=(const LatOpExpr<Op, Lattice<T, D>> &expr)
{
    expr.first.eval(*this, expr.second);

    return *this;
}

// dump to stdout //////////////////////////////////////////////////////////////
template <typename T, unsigned int D>
void Lattice<T, D>::dump(void)
{
    std::string buf;
    Coord<D>    x;

    for (unsigned int i = 0; i < layout_->getLocalVolume(); ++i)
    {
        x   = layout_->getCoord(i);
        buf = "[";
        for (unsigned int d = 0; d < D; ++d)
        {
            buf += strFrom(x[d]) + ((d == D - 1) ? "" : ", ");
        }
        buf += "](" + strFrom(i) + ")= " + strFrom((*this)(i));
        nodeLog(buf);
    }
}

// reduction ///////////////////////////////////////////////////////////////////
template <typename T, unsigned int D>
T Lattice<T, D>::reduce(const T &x, std::function<T(const Vec<T> &buf)> &f)
{
    MPI_Allgather(const_cast<T *>(&x), 1, mpiSiteType_,
                  reduceBuffer_.data(), 1, mpiSiteType_,
                  layout_->getCommGrid());

    return f(reduceBuffer_);
}

template <typename T, unsigned int D>
T Lattice<T, D>::reduceSum(const T &x)
{
    std::function<T(const Vec<T> &buf)> sum;

    sum = [](const Vec<T> &buf)->T{return buf.sum();};

    return reduce(x, sum);
}

// lattice integration /////////////////////////////////////////////////////////
template <typename T, unsigned int D>
T Lattice<T, D>::sum(void)
{
    T localSum;

    localSum = (*this)(0) - (*this)(0);
    FOR_SITE(*this, i)
    {
        localSum += (*this)(i);
    }

    return reduceSum(localSum);
}

template <typename T, unsigned int D>
Vec<T> Lattice<T, D>::sumSlice(const unsigned d)
{
    Vec<T>       localSum(layout_->getLocalDim(d)), sum(layout_->getDim(d));
    T            buf, zero = (*this)(0) - (*this)(0);
    unsigned int xi = layout_->getFirstSite()[d];
    unsigned int xf = xi + layout_->getLocalDim(d);

    localSum.fill(zero);
    for (unsigned int x = 0; x < layout_->getLocalDim(d); ++x)
    {
        FOR_SITE_IN_PLANE(*this, d, x, i)
        {
            localSum(x) += (*this)(i);
        }
    }
    for (unsigned int x = 0; x < layout_->getDim(d); ++x)
    {
        buf    = ((x >= xi)&&(x < xf)) ? localSum(x - xi) : zero;
        sum(x) = reduceSum(buf);
    }

    return sum;
}

END_LATSIM_NAMESPACE

#endif // LatSim_Lattice_hpp_
