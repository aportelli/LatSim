/*
 * LatOp.hpp, part of LatSim
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

#ifndef LatSim_LatOp_hpp_
#define LatSim_LatOp_hpp_

#include <LatSim/Global.hpp>
#include <LatSim/Expressions.hpp>
#include <LatSim/Lattice.hpp>
#include <utility>

BEGIN_NAMESPACE

// forward finite difference ///////////////////////////////////////////////////
template <typename Lat>
class ForwardDiff
{
public:
    explicit ForwardDiff(const unsigned int dir);
    void strong_inline eval(Lat &res, Lat &arg) const;
private:
    unsigned int dir_;
};

template <typename Lat>
ForwardDiff<Lat>::ForwardDiff(const unsigned int dir)
: dir_(dir)
{}

template <typename Lat>
void strong_inline ForwardDiff<Lat>::eval(Lat &res, Lat &arg) const
{
    const unsigned int     localVol = arg.getLayout().getLocalVolume();
    constexpr unsigned int cacheSize = 1024;
    unsigned int           endIndex;
    typename Lat::SiteType cache[cacheSize];

    arg.gatherStart(dir_);
    res = arg;
    arg.gatherWait(dir_);
    for (unsigned int i = 0; i < localVol; i += cacheSize)
    {
        endIndex = (i + cacheSize < localVol) ? cacheSize : (localVol - i - 1);
        for (unsigned int j = 0; j < endIndex; ++j)
        {
            cache[j] = arg(j, dir_);
        }
        for (unsigned int j = 0; j < endIndex; ++j)
        {
            res(j) -= cache[j];
        }
    }
}

template <typename Lat>
LatOpExpr<ForwardDiff<Lat>, Lat> forwardDiff(Lat &l, const unsigned int dir)
{
    return LatOpExpr<ForwardDiff<Lat>, Lat>(ForwardDiff<Lat>(dir), l);
}

// Laplacian ///////////////////////////////////////////////////////////////////
template <typename T, unsigned int D>
class Laplacian
{
public:
    Laplacian(void) = default;
    void strong_inline eval(Lattice<T, D> &res, Lattice<T, D> &arg) const;
};

template <typename T, unsigned int D>
void strong_inline Laplacian<T, D>::eval(Lattice<T, D> &res,
                                         Lattice<T, D> &arg) const
{
    const unsigned int     localVol = arg.getLayout().getLocalVolume();
    unsigned int           endIndex;
    constexpr unsigned int cacheSize = 1024;

    for (unsigned int d = 0; d < 2*D; ++d)
    {
        arg.gatherStart(d);
    }
    for (unsigned int i = 0; i < localVol; ++i)
    {
        res(i) = -2u*D*arg(i);
    }
    for (unsigned int d = 0; d < 2*D; ++d)
    {
        arg.gatherWait(d);
    }
    for (unsigned int d = 0; d < 2*D; ++d)
    {
        T  cache[cacheSize];

        for (unsigned int i = 0; i < localVol; i += cacheSize)
        {
            endIndex = (i + cacheSize < localVol) ? cacheSize
                                                  : (localVol - i - 1);
            for (unsigned int j = 0; j < endIndex; ++j)
            {
                cache[j] = arg(j, d);
            }
            for (unsigned int j = 0; j < endIndex; ++j)
            {
                res(j) += cache[j];
            }
        }
    }
}

template <typename T, unsigned int D>
LatOpExpr<Laplacian<T, D>, Lattice<T, D>> laplacian(Lattice<T, D> &l)
{
    return LatOpExpr<Laplacian<T, D>, Lattice<T, D>>(Laplacian<T, D>(), l);
}

END_NAMESPACE

#endif // LatSim_LatOp_hpp_
