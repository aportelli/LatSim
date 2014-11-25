/*
 * Expressions.hpp, part of LatSim
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

#ifndef LatSim_Expressions_hpp_
#define LatSim_Expressions_hpp_

#include <LatSim/Global.hpp>
#include <tuple>

#define BEGIN_EXPR_NAMESPACE namespace Expr {
#define END_EXPR_NAMESPACE   }

BEGIN_NAMESPACE
BEGIN_EXPR_NAMESPACE

// compile time stack machine, strongly inspired by:
// http://stackoverflow.com/questions/11809052/expression-templates-and-c11

/******************************************************************************
 *                            site tuple factory                              *
 ******************************************************************************/
// site extractor class ////////////////////////////////////////////////////////
template <typename SiteTuple, typename LatTuple, int j, bool end>
class Site
{
public:
    template <typename... Ts>
    static SiteTuple site(int i, const LatTuple& t, Ts && ... args);
};

template <typename SiteTuple, typename LatTuple, int j, bool end>
template <typename... Ts>
SiteTuple
Site<SiteTuple, LatTuple, j, end>::site(int i, const LatTuple& t,
                                        Ts && ... args)
{
    return Site<SiteTuple, LatTuple,  j+1,
    std::tuple_size<LatTuple>::value == j+1>::
    site(i, t , std::forward<Ts>(args)..., std::get<j>(t)[i]);
}

// partial specialization for the last element /////////////////////////////////
template <typename SiteTuple, typename LatTuple, int j>
class Site<SiteTuple, LatTuple, j, true>
{
public:
    template <typename ... Ts>
    static SiteTuple site(int i, LatTuple const& t, Ts && ... args);
};

template <typename SiteTuple, typename LatTuple, int N>
template <typename... Ts>
SiteTuple
Site<SiteTuple, LatTuple, N, true>::site(int i __unused,
                                         const LatTuple& t __unused,
                                         Ts && ... args)
{
    return SiteTuple(std::forward<Ts>(args)...);
}

// function to build a site tuple from a lattice tuple /////////////////////////
template <typename... Lats>
std::tuple<typename std::remove_reference<Lats>::type::SiteType...>
getSites(int i , const std::tuple<Lats...>& latTuple)
{
    return Site<
    std::tuple<typename std::remove_reference<Lats>::type::SiteType...>,
        std::tuple<Lats...>, 0,
        std::tuple_size<std::tuple<Lats...>>::value == 0>::site(i, latTuple);
}

/******************************************************************************
 *                          template stack machine                            *
 ******************************************************************************/
// evaluation of one step in the stack
template <int OpPos, int ArgPos, bool end>
class StackMachine
{
public:
    template <typename... Ts, typename... Ops>
    static void eval(std::tuple<Ts...> &args, const std::tuple<Ops...> &ops);
};

template <int OpPos, int ArgPos, bool end>
template <typename... Ts, typename... Ops>
void StackMachine<OpPos, ArgPos, end>::eval(std::tuple<Ts...> &args,
                                            const std::tuple<Ops...> &ops)
{
    auto           &currArg = std::get<ArgPos>(args);
    auto           &nextArg = std::get<ArgPos + 1>(args);
    auto           &currOp  = std::get<OpPos>(ops);
    constexpr bool isLast   = (sizeof...(Ops) == OpPos + 1);

    nextArg = currOp.eval(currArg, nextArg);
    StackMachine<OpPos + 1, ArgPos + 1, isLast>::eval(args, ops);
}

// specialization for the stack's end
template <int OpPos, int ArgPos>
class StackMachine<OpPos, ArgPos, true>
{
public:
    template <typename... Ts, typename... Ops>
    static void eval(std::tuple<Ts...>& args __unused,
                     const std::tuple<Ops...> & ops __unused)
    {}
};

// global evaluation function
template <typename T, typename... Ts, typename... Ops>
T eval(const std::tuple<Ts...>& args, const std::tuple<Ops...> & ops)
{
    StackMachine<0, 0, false>::eval(const_cast<std::tuple<Ts...> &>(args), ops);

    return std::get<sizeof...(Ops)>(args);
}

/******************************************************************************
 *                          elementary operations                             *
 ******************************************************************************/
template <typename LhsT, typename RhsT>
class Add
{
public:
    static auto eval(const LhsT &lhs, const RhsT &rhs)->decltype(lhs + rhs)
    {
        return lhs + rhs;
    }
};

template <typename LhsT, typename RhsT>
class RSub
{
public:
    static auto eval(const LhsT &lhs, const RhsT &rhs)->decltype(lhs - rhs)
    {
        return lhs - rhs;
    }
};

template <typename LhsT, typename RhsT>
class LSub
{
public:
    static auto eval(const LhsT &lhs, const RhsT &rhs)->decltype(lhs - rhs)
    {
        return rhs - lhs;
    }
};

template <typename LhsT, typename RhsT>
class RMul
{
public:
    static auto eval(const LhsT &lhs, const RhsT &rhs)->decltype(lhs*rhs)
    {
        return lhs*rhs;
    }
};

template <typename LhsT, typename RhsT>
class LMul
{
public:
    static auto eval(const LhsT &lhs, const RhsT &rhs)->decltype(lhs*rhs)
    {
        return rhs*lhs;
    }
};

template <typename LhsT, typename RhsT>
class RDiv
{
public:
    static auto eval(const LhsT &lhs, const RhsT &rhs)->decltype(lhs/rhs)
    {
        return lhs/rhs;
    }
};

template <typename LhsT, typename RhsT>
class LDiv
{
public:
    static auto eval(const LhsT &lhs, const RhsT &rhs)->decltype(lhs/rhs)
    {
        return rhs/lhs;
    }
};

END_EXPR_NAMESPACE
END_NAMESPACE

#endif // LatSim_Expressions_hpp_
