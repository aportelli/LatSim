/*
 * Phi4.hpp, part of LatSim
 *
 * Copyright (C) 2015 Antonin Portelli
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

#ifndef LatSim_Phi4_hpp_
#define LatSim_Phi4_hpp_

#include <LatSim/Global.hpp>
#include <LatSim/Lattice.hpp>
#include <random>

#define BEGIN_PHI4_NAMESPACE namespace Phi4 {
#define END_PHI4_NAMESPACE   }

BEGIN_LATSIM_NAMESPACE

BEGIN_PHI4_NAMESPACE

template <unsigned int D>
inline void momGen(Lattice<double, D> &phi __unused, Lattice<double, D> &pi)
{
    std::normal_distribution<double> d;
    const Layout<D>                  &layout = pi.getLayout();

    FOR_SITE(pi, i)
    {
        pi(i) = d(layout.getRng(i));
    }
}

template <unsigned int D>
inline void moveField(Lattice<double, D> &phi, Lattice<double, D> &pi,
                      const double step)
{
    FOR_SITE(phi, i)
    {
        phi(i) += pi(i)*step;
    }
}

struct ActionParam
{
    double lambda, kappa;
};

template <unsigned int D>
class Hamiltonian
{
public:
    // constructor
    Hamiltonian(const ActionParam &param): param_(param) {};
    // destructor
    virtual ~Hamiltonian(void) = default;
    // evaluation
    inline double operator()(Lattice<double, D> &phi,
                             Lattice<double, D> &pi) const
    {
        const double lambda = param_.lambda, kappa = param_.kappa;
        double       h = 0., dSum, phi2;

        for (unsigned int d = 0; d < D; ++d)
        {
            phi.gather(d);
        }
        FOR_SITE(phi, i)
        {
            dSum = 0.;
            for (unsigned int d = 0; d < D; ++d)
            {
                dSum += phi(i, d);
            }
            phi2  = phi(i)*phi(i);
            h    += -2*kappa*dSum*phi(i) + phi2 + lambda*(phi2-1.)*(phi2-1.);
            h    += .5*pi(i)*pi(i);
        }

        return h;
    }
private:
    const ActionParam &param_;
};

template <unsigned int D>
class MoveMom
{
public:
    // constructor
    MoveMom(const ActionParam &param): param_(param) {};
    // destructor
    virtual ~MoveMom(void) = default;
    // evaluation
    inline void operator()(Lattice<double, D> &phi, Lattice<double, D> &pi,
                           const double step) const
    {
        const double lambda = param_.lambda, kappa = param_.kappa;
        double       dSum, force;

        for (unsigned int d = 0; d < 2*D; ++d)
        {
            phi.gather(d);
        }
        FOR_SITE(phi, i)
        {
            dSum = 0.;
            for (unsigned int d = 0; d < 2*D; ++d)
            {
                dSum += phi(i, d);
            }
            force  = 2.*kappa*dSum - 2.*phi(i)
                     - lambda*4.*(phi(i)*phi(i)-1.)*phi(i);
            pi(i) += force*step;
        }
    }
private:
    const ActionParam &param_;
};

END_PHI4_NAMESPACE

END_LATSIM_NAMESPACE

#endif // LatSim_Phi4_hpp_
