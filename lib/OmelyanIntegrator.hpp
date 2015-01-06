/*
 * OmelyanIntegrator.hpp, part of LatSim
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

#ifndef LatSim_OmelyanIntegrator_hpp_
#define LatSim_OmelyanIntegrator_hpp_

#include <LatSim/Global.hpp>
#include <LatSim/Integrator.hpp>

BEGIN_NAMESPACE

/******************************************************************************
 *                           OmelyanIntegrator                                *
 ******************************************************************************/
template <typename... Fs>
class OmelyanIntegrator: public Integrator<Fs...>
{
public:
    typedef typename Integrator<Fs...>::EvolFunc EvolFunc;
public:
    static constexpr double xi = 0.1931833;
public:
    // constructor
    OmelyanIntegrator(const EvolFunc &moveFields, const EvolFunc &moveMoms);
    // destructor
    virtual ~OmelyanIntegrator(void) = default;
    // molecular dynamic evolution
    virtual void evolve(Fs &... fields, const double step,
                        const unsigned int nStep);
private:

};

/******************************************************************************
 *                      OmelyanIntegrator implementation                      *
 ******************************************************************************/
template <typename... Fs>
OmelyanIntegrator<Fs...>::OmelyanIntegrator(const EvolFunc &moveFields,
                                            const EvolFunc &moveMoms)
: Integrator<Fs...>(moveFields, moveMoms)
{}

template <typename... Fs>
void OmelyanIntegrator<Fs...>::evolve(Fs &... fields, const double step,
                                      const unsigned int nStep)
{
    for (unsigned int i = 0; i < nStep; ++i)
    {
        this->moveFields(fields..., step*xi);
        this->moveMoms(fields..., step/2.);
        this->moveFields(fields..., step*(1.-2.*xi));
        this->moveMoms(fields..., step/2.);
        this->moveFields(fields..., step*xi);
    }
}

END_NAMESPACE

#endif // LatSim_OmelyanIntegrator_hpp_
