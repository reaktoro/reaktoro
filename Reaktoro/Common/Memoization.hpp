// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// C++ includes
#include <map>
#include <memory>
#include <tuple>

// Reaktoro includes
#include <Reaktoro/Common/TraitsUtils.hpp>

namespace Reaktoro {

template<typename Ret, typename... Args>
auto memoize(std::function<Ret(Args...)> f) -> std::function<Ret(Args...)>
{
    auto cache = std::make_shared<std::map<std::tuple<Args...>, Ret>>();
    return [=](Args... args) mutable -> Ret
    {
        std::tuple<Args...> t(args...);
        if(cache->find(t) == cache->end())
            (*cache)[t] = f(args...);
        return (*cache)[t];
    };
}

template<typename Fun, EnableIf<!isFunction<Fun>>...>
auto memoize(Fun f)
{
    return memoize(std::function(f));
}

template<typename Ret, typename... Args>
auto memoizeLast(std::function<Ret(Args...)> f) -> std::function<Ret(Args...)>
{
    std::tuple<typename std::decay<Args>::type...> cache;
    Ret result = Ret();
    auto firsttime = true;
    return [=](Args... args) mutable -> Ret
    {
        if(std::tie(args...) == cache && !firsttime)
            return Ret(result);
        cache = std::make_tuple(args...);
        firsttime = false;
        return result = f(args...);
    };
}

template<typename Fun, EnableIf<!isFunction<Fun>>...>
auto memoizeLast(Fun f)
{
    return memoizeLast(std::function(f));
}

template<typename Ret, typename... Args>
auto memoizeLastPtr(std::function<Ret(Args...)> f) -> std::shared_ptr<std::function<Ret(Args...)>>
{
    return std::make_shared<std::function<Ret(Args...)>>(memoizeLast(f));
}

template<typename Ret, typename... Args>
auto dereference(const std::shared_ptr<std::function<Ret(Args...)>>& f) -> std::function<Ret(Args...)>
{
    return [=](Args... args) -> Ret { return (*f)(args...); };
}

} // namespace Reaktoro

