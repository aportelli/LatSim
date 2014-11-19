/*
 * Logger.hpp, part of LatSim
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

#ifndef LatSim_Logger_hpp_
#define LatSim_Logger_hpp_

#include <LatSim/Global.hpp>

BEGIN_NAMESPACE

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
    void nodeLog(const std::string &msg) const;
    void nodeLog(const std::string &&msg) const;
private:
    // check
    void checkMpi(void) const;
private:
    std::string name_;
};

END_NAMESPACE

#endif // LatSim_Logger_hpp_
