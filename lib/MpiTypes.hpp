/*
 * MpiTypes.hpp, part of LatSim
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

#ifndef LatSim_MpiTypes_hpp_
#define LatSim_MpiTypes_hpp_

#include <LatSim/Global.hpp>

BEGIN_NAMESPACE

// MPI type constructor is embeded inside a class to allow partial template
// specialization

// generic case
template <typename T>
class MpiType
{
public:
    static void make(MPI_Datatype &type __unused)
    {
        locGlobalError("MPI datatype not implemented");
    }
};

// unsigned int
template <>
class MpiType<unsigned int>
{
public:
    static void make(MPI_Datatype &type)
    {
        MPI_Type_contiguous(1, MPI_UNSIGNED, &type);
    }
};

// float
template <>
class MpiType<float>
{
public:
    static void make(MPI_Datatype &type)
    {
        MPI_Type_contiguous(1, MPI_FLOAT, &type);
    }
};

// double
template <>
class MpiType<double>
{
public:
    static void make(MPI_Datatype &type)
    {
        MPI_Type_contiguous(1, MPI_DOUBLE, &type);
    }
};

// matrix of float
template <int nRow, int nCol>
class MpiType<SFMat<nRow, nCol>>
{
public:
    static void make(MPI_Datatype &type)
    {
        MPI_Type_contiguous(nRow*nCol, MPI_FLOAT, &type);
    }
};

// matrix of double
template <int nRow, int nCol>
class MpiType<SDMat<nRow, nCol>>
{
public:
    static void make(MPI_Datatype &type)
    {
        MPI_Type_contiguous(nRow*nCol, MPI_DOUBLE, &type);
    }
};

END_NAMESPACE

#endif // LatSim_MpiTypes_hpp_
