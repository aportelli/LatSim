/*
 * VerletIntegrator.hpp, part of LatSim
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

#ifndef LatSim_VerletIntegrator_hpp_
#define LatSim_VerletIntegrator_hpp_

#include <LatSim/Global.hpp>
#include <LatSim/Integrator.hpp>

BEGIN_NAMESPACE

/******************************************************************************
 *                           VerletIntegrator                                 *
 ******************************************************************************/
template <typename... Fs>
class VerletIntegrator: public Integrator<Fs...>
{
public:
    typedef typename Integrator<Fs...>::EvolFunc EvolFunc;
public:
    // constructor
    VerletIntegrator(const EvolFunc &moveFields, const EvolFunc &moveMoms);
    // destructor
    virtual ~VerletIntegrator(void) = default;
    // molecular dynamic evolution
    virtual void evolve(Fs &... fields, const double step,
                        const unsigned int nStep);
};

/******************************************************************************
 *                      VerletIntegrator implementation                       *
 ******************************************************************************/
template <typename... Fs>
VerletIntegrator<Fs...>::VerletIntegrator(const EvolFunc &moveFields,
                                          const EvolFunc &moveMoms)
: Integrator<Fs...>(moveFields, moveMoms)
{}

template <typename... Fs>
void VerletIntegrator<Fs...>::evolve(Fs &... fields, const double step,
                                     const unsigned int nStep)
{
    for (unsigned int i = 0; i < nStep; ++i)
    {
        this->moveFields(fields..., step/2.);
        this->moveMoms(fields..., step);
        this->moveFields(fields..., step/2.);
    }
}

END_NAMESPACE

#endif // LatSim_VerletIntegrator_hpp_
