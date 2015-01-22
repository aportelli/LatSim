/*
 * Global.hpp, part of LatSim
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

#ifndef LatSim_Global_hpp_
#define	LatSim_Global_hpp_

#include <array>
#include <complex>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <sstream>
#include <type_traits>
#include <vector>
#include <cstdlib>

//#define EIGEN_NO_DEBUG
#pragma GCC diagnostic ignored "-Wdeprecated-register"
#pragma GCC diagnostic ignored "-Wconversion"

#define EIGEN_DONT_PARALLELIZE
#include <LatSim/Eigen/Dense>

#pragma GCC diagnostic warning "-Wdeprecated-register"
#pragma GCC diagnostic warning "-Wconversion"

#include <mpi.h>

#ifdef LATSIM_FORCE_INLINE
#define strong_inline __attribute__((always_inline)) inline
#else
#define strong_inline inline
#endif

#ifdef LATSIM_FLATTEN
#define flatten __attribute__((flatten))
#else
#define flatten
#endif

#define BEGIN_NAMESPACE namespace LatSim {
#define END_NAMESPACE }

// macro utilities
#define unique_arg(...) __VA_ARGS__
#define DEBUG_VAR(x) std::cout << #x << "= "  << x << std::endl
#define DEBUG_MAT(m) std::cout << #m << "=\n" << m << std::endl

// attribute to switch off unused warnings with gcc
#if (!defined(__GNUC__))||\
    (defined(__clang__)&&(__clang_major__ <= 3)&&(__clang_minor__ <= 4))
#define __unused
#endif

// copy/assignement from Eigen expression
#define EIGEN_EXPR_CTOR(ctorName, Class, Base, ExprType) \
template <typename Derived>\
ctorName(const ExprType<Derived> &m): Base(m) {}\
template<typename Derived>\
Class & operator=(const ExprType<Derived> &m)\
{\
    this->Base::operator=(m);\
    return *this;\
}

BEGIN_NAMESPACE

// lattice coordinate type /////////////////////////////////////////////////////
template <unsigned int D>
using Coord = std::array<unsigned int, D>;

// lattice base type ///////////////////////////////////////////////////////////
class LatticeObj {};

// Eigen type aliases //////////////////////////////////////////////////////////
const int dynamic = -1;

// array types
template <typename Derived>
using ArrayExpr = Eigen::ArrayBase<Derived>;

template <typename T, int nRow = dynamic, int nCol = dynamic>
using Array = Eigen::Array<T, nRow, nCol>;

// matrix types
template <typename Derived>
using MatExpr = Eigen::MatrixBase<Derived>;

template <typename T, int nRow = dynamic, int nCol = dynamic>
using Mat = Eigen::Matrix<T, nRow, nCol>;

template <int nRow, int nCol>
using SFMat = Eigen::Matrix<float, nRow, nCol>;

template <int nRow, int nCol>
using SDMat = Eigen::Matrix<double, nRow, nCol>;

template <int nRow, int nCol>
using SCMat = Eigen::Matrix<std::complex<double>, nRow, nCol>;

typedef SFMat<dynamic, dynamic> FMat;
typedef SDMat<dynamic, dynamic> DMat;
typedef SCMat<dynamic, dynamic> CMat;

// vector types
template <typename T, int size = dynamic>
using Vec = Mat<T, size, 1>;

template <int size>
using SIVec = Vec<int, size>;

template <int size>
using SUVec = Vec<unsigned int, size>;

template <int size>
using SFVec = Vec<float, size>;

template <int size>
using SDVec = Vec<double, size>;

template <int size>
using SCVec = Vec<std::complex<double>, size>;

typedef SIVec<dynamic> IVec;
typedef SUVec<dynamic> UVec;
typedef SDVec<dynamic> DVec;
typedef SCVec<dynamic> CVec;

#define FOR_VEC(vec, i)  for (LatSim::Index i = 0; i < (vec).size(); ++i)
#define FOR_ARRAY(ar, i) FOR_VEC(ar, i)

// block types
template <typename Derived>
using Block      = Eigen::Block<Derived>;
template <typename Derived>
using ConstBlock = const Eigen::Block<const Derived>;

template <typename Derived>
using Row      = typename Derived::RowXpr;
template <typename Derived>
using ConstRow = typename Derived::ConstRowXpr;

template <typename Derived>
using Col = typename Derived::ColXpr;
template <typename Derived>
using ConstCol = typename Derived::ConstColXpr;

// map type
template <typename Derived>
using Map = Eigen::Map<Derived>;
template <typename Derived>
using ConstMap = Eigen::Map<const Derived>;

// Index type //////////////////////////////////////////////////////////////////
typedef Mat<double>::Index MatIndex;

// Error handling //////////////////////////////////////////////////////////////
#define SRC_LOC LatSim::strFrom(__FUNCTION__) + " at " +\
                LatSim::strFrom(__FILE__) + ":" + LatSim::strFrom(__LINE__)
#define locGlobalError(msg) LatSim::globalError(msg, SRC_LOC)

void globalError(const std::string msg, const std::string loc = "");

// Indexing helpers ////////////////////////////////////////////////////////////
template <size_t D>
inline unsigned int coordToRowMajor(const Coord<D> &x, const Coord<D> &dim)
{
    unsigned int ind;

    ind = dim[1]*x[0];
    for (unsigned int d = 1; d < D - 1; ++d)
    {
        ind = dim[d+1]*(x[d] + ind);
    }
    ind += x[D-1];

    return ind;
}

template <size_t D>
inline Coord<D> rowMajorToCoord(const unsigned int ind, const Coord<D> &dim)
{
    Coord<D>     x;
    unsigned int j, dimprod, ud;

    j       = ind;
    dimprod = 1;
    for (int d = D - 1; d >= 0; --d)
    {
        ud       = static_cast<unsigned int>(d);
        x[ud]    = (j/dimprod)%dim[ud];
        j       -= dimprod*x[ud];
        dimprod *= dim[ud];
    }

    return x;
}

// Environment /////////////////////////////////////////////////////////////////
namespace Env
{
    extern const std::string fullName;
    extern const std::string name;
    extern const std::string version;
    extern const std::string msgPrefix;
    // empty function for library test
    void function(void);
}

// String conversions //////////////////////////////////////////////////////////
template <typename T>
inline T strTo(const std::string &str)
{
    T buf;
    std::istringstream stream(str);
    
    stream >> buf;
    
    return buf;
}

// optimized specializations
template <>
inline float strTo<float>(const std::string &str)
{
    return strtof(str.c_str(), (char **)NULL);
}
template <>
inline double strTo<double>(const std::string &str)
{
    return strtod(str.c_str(), (char **)NULL);
}
template <>
inline int strTo<int>(const std::string &str)
{
    return (int)(strtol(str.c_str(), (char **)NULL, 10));
}
template <>
inline long strTo<long>(const std::string &str)
{
    return strtol(str.c_str(), (char **)NULL, 10);
}
template <>
inline std::string strTo<std::string>(const std::string &str)
{
    return str;
}

template <typename T>
inline std::string strFrom(const T x)
{
    std::ostringstream stream;
    
    stream << x;
    
    return stream.str();
}

// specialization for vectors
template<>
inline DVec strTo<DVec>(const std::string &str)
{
    DVec                res;
    std::vector<double> vbuf;
    double              buf;
    std::istringstream  stream(str);
    
    while (!stream.eof())
    {
        stream >> buf;
        vbuf.push_back(buf);
    }
    res = Map<DVec>(vbuf.data(), static_cast<MatIndex>(vbuf.size()));
    
    return res;
}

template<>
inline IVec strTo<IVec>(const std::string &str)
{
    IVec                res;
    std::vector<int>    vbuf;
    int                 buf;
    std::istringstream  stream(str);
    
    while (!stream.eof())
    {
        stream >> buf;
        vbuf.push_back(buf);
    }
    res = Map<IVec>(vbuf.data(), static_cast<MatIndex>(vbuf.size()));
    
    return res;
}

template<>
inline UVec strTo<UVec>(const std::string &str)
{
    UVec                      res;
    std::vector<unsigned int> vbuf;
    unsigned int              buf;
    std::istringstream        stream(str);

    while (!stream.eof())
    {
        stream >> buf;
        vbuf.push_back(buf);
    }
    res = Map<UVec>(vbuf.data(), static_cast<MatIndex>(vbuf.size()));

    return res;
}

END_NAMESPACE

#endif // LatSim_Global_hpp_
