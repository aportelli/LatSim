/*
 * Integrator.hpp, part of LatSim
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

#ifndef LatSim_Integrator_hpp_
#define LatSim_Integrator_hpp_

#include <LatSim/Global.hpp>
#include <functional>

BEGIN_NAMESPACE

/******************************************************************************
 *                               Integrator                                   *
 ******************************************************************************/
template <typename... Fs>
class Integrator
{
public:
    typedef std::function<void(Fs &..., const double)> EvolFunc;
public:
    // constructor
    Integrator(const EvolFunc &moveFields, const EvolFunc &moveMoms);
    // destructor
    virtual ~Integrator(void) = default;
    // integrator steps
    void moveFields(Fs &... fields, const double step) const;
    void moveMoms(Fs &... fields, const double step) const;
    // molecular dynamic evolution
    virtual void evolve(Fs &... fields, const double step,
                        const unsigned int nStep) const = 0;
private:
    EvolFunc moveFields_, moveMoms_;
};

/******************************************************************************
 *                       IntegratorBase implementation                        *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename... Fs>
Integrator<Fs...>::Integrator(const EvolFunc &moveFields,
                              const EvolFunc &moveMoms)
: moveFields_(moveFields)
, moveMoms_(moveMoms)
{}

// integrator steps ////////////////////////////////////////////////////////////
template <typename... Fs>
void Integrator<Fs...>::moveFields(Fs &... fields, const double step) const
{
    moveFields_(fields..., step);
}

template <typename... Fs>
void Integrator<Fs...>::moveMoms(Fs &... fields, const double step) const
{
    moveMoms_(fields..., step);
}

END_NAMESPACE

#endif // LatSim_Integrator_hpp_
