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
#include <LatSim/Layout.hpp>

BEGIN_NAMESPACE

/******************************************************************************
 *                                Lattice                                     *
 ******************************************************************************/

template <typename T, unsigned long D>
class Lattice
{
public:
    // constructor
    explicit Lattice(const LayoutObject *layout = globalLayout);
    // destructor
    virtual ~Lattice(void) = default;
private:
    const Layout<D>                   *layout_;
    std::unique_ptr<T>                data_;
    std::array<std::unique_ptr<T>, D> commBuffer_;
};

/******************************************************************************
 *                          Lattice implementation                            *
 ******************************************************************************/
template <typename T, unsigned long D>
Lattice<T, D>::Lattice(const LayoutObject *layout)
{
    layout_ = dynamic_cast<const Layout<D> *>(layout);
    data_.reset(new T[layout_->getLocalVolume()]);
    for (unsigned int d = 0; d < D; ++d)
    {
        commBuffer_[d].reset(new T[layout_->getLocalSurface(d)]);
    }
}

END_NAMESPACE

#endif // LatSim_Lattice_hpp_
