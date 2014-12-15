/*
 * Global.cpp, part of LatSim
 *
 * Copyright (C) 2013 - 2014 Antonin Portelli
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

#include <LatSim/Global.hpp>
#include <LatSim/includes.hpp>

using namespace std;
using namespace LatSim;

const string Env::fullName  = PACKAGE_STRING;
const string Env::name      = PACKAGE_NAME;
const string Env::version   = PACKAGE_VERSION;
const string Env::msgPrefix = "[" + strFrom(PACKAGE_NAME) + " v"
                              + strFrom(PACKAGE_VERSION) + "] ";

void LatSim::globalError(const string msg, const string loc)
{
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        cerr << "GLOBAL ERROR: " << msg;
        if (!loc.empty())
        {
            cerr << " (" << loc << ")" << endl;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
    {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
}

void Env::function(void)
{}
