/*
 * Logger.hpp, part of LatSim
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

#ifndef LatSim_Logger_hpp_
#define LatSim_Logger_hpp_

#include <LatSim/Global.hpp>
#include <iomanip>
#include <cstdarg>

BEGIN_LATSIM_NAMESPACE

/******************************************************************************
 *                                  Logger                                    *
 ******************************************************************************/
class Logger
{
public:
    static const unsigned int width = 20;
public:
    // constructor
    explicit Logger(const std::string name);
    // destructor
    virtual ~Logger(void) = default;
    // IO
    void masterLog(const std::string &msg) const;
    void masterLog(const std::string &&msg) const;
    void masterLog(const char *fmt, ...) const;
    void nodeLog(const std::string &msg) const;
    void nodeLog(const std::string &&msg) const;
private:
    // check
    void checkMpi(void) const;
private:
    std::string name_;
};

/******************************************************************************
 *                           Logger implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
Logger::Logger(const std::string name)
: name_(name)
{}

// IO //////////////////////////////////////////////////////////////////////////
void Logger::masterLog(const std::string &msg) const
{
    int rank;

    checkMpi();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        std::cout << std::left << std::setw(width) << "[" + name_ + "]";
        std::cout << ": " << msg << std::endl;
    }
}

void Logger::masterLog(const std::string &&msg) const
{
    masterLog(msg);
}

void Logger::masterLog(const char *fmt, ...) const
{
    static char  buf[128];
    va_list args;

    va_start(args, fmt);
    vsprintf(buf, fmt, args);
    va_end(args);
    masterLog(std::string(buf));
}

void Logger::nodeLog(const std::string &msg) const
{
    int rank, size;

    checkMpi();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    for (int r = 0; r < size; ++r)
    {
        if (r == rank)
        {
            std::cout << std::left << std::setw(width);
            std::cout << "{" + name_ + "(" + strFrom(rank) + ")}" << ": ";
            std::cout << msg << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void Logger::nodeLog(const std::string &&msg) const
{
    nodeLog(msg);
}

// check ///////////////////////////////////////////////////////////////////////
void Logger::checkMpi(void) const
{
    int isInit;

    MPI_Initialized(&isInit);
    if (!isInit)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        locGlobalError("communications not initialized, please declare"
                       " a Layout object");
    }
}

END_LATSIM_NAMESPACE

#endif // LatSim_Logger_hpp_
