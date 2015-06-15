/*
 * Expressions.hpp, part of LatSim
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

#ifndef LatSim_Expressions_hpp_
#define LatSim_Expressions_hpp_

#include <LatSim/Global.hpp>
#include <tuple>

#define BEGIN_EXPR_NAMESPACE namespace Expr {
#define END_EXPR_NAMESPACE   }

BEGIN_LATSIM_NAMESPACE

template <typename Op, typename... Ts>
class LatExpr: public std::tuple<Op, Ts...>, public LatticeObj
{
public:
    using std::tuple<Op, Ts...>::tuple;
};

template <typename Op, typename Lat>
using LatOpExpr = std::pair<Op, Lat &>;

BEGIN_EXPR_NAMESPACE

template<int... is>
class ISeq {};

template<int N, int... is>
class SeqGen: public SeqGen<N-1, N-1, is...> {};

template<int... is>
class SeqGen<1, is...>
{
public:
    static strong_inline ISeq<is...> seq(void) {return ISeq<is...>();};
};

template <typename T>
auto strong_inline eval(const unsigned int i, const T &arg)->decltype(arg(i))
{
    return arg(i);
}

template <typename Op, typename... Ts, int... is>
auto strong_inline eval(const unsigned int i,
                        const LatExpr<Op, Ts...> &expr,
                        const ISeq<is...> &seq __dumb)
->decltype(std::get<0>(expr).eval(eval(i, std::get<is>(expr))...))
{
    return std::get<0>(expr).eval(eval(i, std::get<is>(expr))...);
}

template <typename Op, typename... Ts>
auto strong_inline eval(const unsigned int i,
                        const LatExpr<Op, Ts...> &expr)
->decltype(eval(i, expr, SeqGen<sizeof...(Ts)+1>::seq()))
{
    return eval(i, expr, SeqGen<sizeof...(Ts)+1>::seq());
}

/******************************************************************************
 *                          elementary operations                             *
 ******************************************************************************/
template <typename LhsT, typename RhsT>
class Add
{
public:
    static auto strong_inline eval(const LhsT &lhs, const RhsT &rhs)
    ->decltype(lhs + rhs)
    {
        return lhs + rhs;
    }
};

template <typename LhsT, typename RhsT>
class Sub
{
public:
    static auto strong_inline eval(const LhsT &lhs, const RhsT &rhs)
    ->decltype(lhs - rhs)
    {
        return lhs - rhs;
    }
};

template <typename LhsT, typename RhsT>
class Mul
{
public:
    static auto strong_inline eval(const LhsT &lhs, const RhsT &rhs)
    ->decltype(lhs*rhs)
    {
        return lhs*rhs;
    }
};

template <typename LhsT, typename RhsT>
class Div
{
public:
    static auto strong_inline eval(const LhsT &lhs, const RhsT &rhs)
    ->decltype(lhs/rhs)
    {
        return lhs/rhs;
    }
};

END_EXPR_NAMESPACE
END_LATSIM_NAMESPACE

#endif // LatSim_Expressions_hpp_
