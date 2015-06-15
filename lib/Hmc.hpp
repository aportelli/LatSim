/*
 * Hmc.hpp, part of LatSim
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

#ifndef LatSim_Hmc_hpp_
#define LatSim_Hmc_hpp_

#include <LatSim/Global.hpp>
#include <LatSim/Integrator.hpp>
#include <functional>
#include <random>
#include <utility>

BEGIN_LATSIM_NAMESPACE

/******************************************************************************
 *                                  Hmc                                       *
 ******************************************************************************/
template <typename... Fs>
class Hmc: public Logger
{
public:
    typedef std::function<double(Fs &...)> Functional;
    typedef std::function<void(Fs &...)>   MomGenFunc;
    typedef Integrator<Fs...>              IntegratorType;

public:
    // constructor
    Hmc(const Functional &hamiltonian, const MomGenFunc &momGen,
        const IntegratorType &integrator, const unsigned trajNStep,
        const unsigned int trajPretherm = 50, const double trajLength = 1.);
    // destructor
    virtual ~Hmc(void) = default;
    // access
    int getTraj(void);
    // HMC step
    bool update(Fs &... fields);
    // HMC statistics
    double getExpDH(void) const;
    double getAcceptRate(void) const;
    void   printStat(void) const;
private:
    // accept-reject
    bool accept(const double dH);
private:
    Functional                             hamiltonian_;
    MomGenFunc                             momGen_;
    const IntegratorType                   &integrator_;
    unsigned int                           trajCount_{0}, trajNStep_;
    unsigned int                           tryCount_{0}, trajPreTherm_;
    double                                 trajStep_, trajLength_, expDH_{0.};
    std::tuple<Fs...>                      initFields_;
    std::mt19937                           gen_;
    std::uniform_real_distribution<double> lin_;
};

/******************************************************************************
 *                            Hmc implementation                              *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename... Fs>
Hmc<Fs...>::Hmc(const Functional &hamiltonian, const MomGenFunc &momGen,
                const IntegratorType &integrator, const unsigned trajNStep,
                const unsigned int trajPretherm, const double trajLength)
: Logger("HMC")
, hamiltonian_(hamiltonian)
, momGen_(momGen)
, integrator_(integrator)
, trajNStep_(trajNStep)
, trajPreTherm_(trajPretherm)
, trajStep_(trajLength/static_cast<double>(trajNStep))
, trajLength_(trajLength)
, lin_(0., 1.)
{
    std::random_device dev;

    gen_.seed(static_cast<std::mt19937::result_type>(dev()));
    masterLog("trajectory length= " + strFrom(trajLength_));
    masterLog("trajectory  nStep= " + strFrom(trajNStep_));
    masterLog("trajectory   step= " + strFrom(trajStep_));
}

// access //////////////////////////////////////////////////////////////////////
template <typename... Fs>
int Hmc<Fs...>::getTraj(void)
{
    return trajCount_ - 1;
}

// HMC step ////////////////////////////////////////////////////////////////////
template <typename... Fs>
bool Hmc<Fs...>::update(Fs &... fields)
{
    double      dH, localDH;
    std::string msg;

    // save initial fields
    initFields_ = std::tuple<Fs &...>(fields...);

    // generate momenta
    momGen_(fields...);

    // compute initial energy
    localDH = hamiltonian_(fields...);

    // evolve
    integrator_.evolve(fields..., trajStep_, trajNStep_);

    // accept/reject
    localDH = hamiltonian_(fields...) - localDH;
    MPI_Allreduce(&localDH, &dH, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    masterLog("exp(-DH)= " + strFrom(exp(-dH)));
    msg  = "attempt " + strFrom(tryCount_);
    msg += (trajCount_ < trajPreTherm_) ? " (pre-thermalization)" : "";
    if (accept(dH))
    {
        msg += ": trajectory " + strFrom(trajCount_) + " accepted";
        masterLog(msg);
        tryCount_++;
        trajCount_++;
        if (trajCount_ == trajPreTherm_)
        {
            masterLog("pre-thermalization done");
            tryCount_ = 0;
        }
        else if (trajCount_ > trajPreTherm_)
        {
            expDH_ += exp(-dH);
        }
        
        return true;
    }
    else
    {
        tryCount_++;
        msg += ": trajectory " + strFrom(trajCount_) + " rejected";
        masterLog(msg);
        std::tie(fields...) = initFields_;

        return false;
    }
}

// HMC statistics //////////////////////////////////////////////////////////////
template <typename... Fs>
double Hmc<Fs...>::getExpDH(void) const
{
    return expDH_/static_cast<double>(trajCount_-trajPreTherm_);
}

template <typename... Fs>
double Hmc<Fs...>::getAcceptRate(void) const
{
    return static_cast<double>(trajCount_-trajPreTherm_)/static_cast<double>(tryCount_);
}

template <typename... Fs>
void Hmc<Fs...>::printStat(void) const
{
    if (trajCount_ > trajPreTherm_)
    {
        masterLog("<exp(-DH)>     = " + strFrom(getExpDH()));
        masterLog("acceptance rate= " + strFrom(getAcceptRate()));
    }
}

// accept-reject ///////////////////////////////////////////////////////////////
template <typename... Fs>
bool Hmc<Fs...>::accept(const double dH)
{
    if (dH < 0)
    {
        return true;
    }
    else
    {
        int    rank;
        double r;

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0)
        {
            r = lin_(gen_);
        }
        MPI_Bcast(&r, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        return (r < exp(-dH));
    }
}

END_LATSIM_NAMESPACE

#endif // LatSim_Hmc_hpp_
