// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

// Reaktoro includes
#include <Reaktoro/Common/Meta.hpp>
#include <Reaktoro/Common/TraitsUtils.hpp>
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {
namespace detail {

/// Type trait used to check if @p a and @p b of type @p T have the same value.
/// Define a template specialization for a different type `U` if `operator==`
/// is not adequate. For example, `operator==` for `Eigen::ArrayBase<Derived>`
/// objects return an expression that needs to be evaluated with methods
/// `all()` or `any()`. Thus, in order for memoization of functions with this
/// array type to compile, a template specialization for
/// `SameValue<Eigen::ArrayBase<Derived>>` is needed.
template<typename T>
struct SameValue
{
    static constexpr auto check(const T& a, const T& b) { return a == b; }
};

/// Type trait used to assing the value of @p b into @p a.
/// Define a template specialization for a different type `U` if `operator=` is
/// not adequate. For example, type `Param` has an underlying shared pointer,
/// with `Param::operator=(Param)` defined in a way that the underlying
/// pointers are assigned instead. However, what we need for memoization is
/// that the value inside the pointers are updated. Otherwise it is not
/// possible to determine if the last arguments in the memoized function call
/// (containing `Param` objects) have changed values or not.
template<typename T>
struct AssignValue
{
    static constexpr auto apply(T& a, const T& b) { a = b; }
};

/// Type trait used to create a clone of value @p a of type @p T.
/// Define a template specialization for a different type `U` if the default
/// behavior is not adequate. For example, type `Param` has an underlying
/// shared pointer and the default clone behavior below will simply perform a
/// shallow copy of this pointer. Instead, the Param::clone method should be
/// used and to achieve this, the specialization `CloneValue<Param>` is needed.
template<typename T>
struct CloneValue
{
    static constexpr auto apply(const T& a) -> T { return a; }
};

/// Return true if @p a and @p b have the same value in the sense of type trait `SameValue<T>` for type @p T.
template<typename T>
constexpr auto sameValue(const T& a, const T& b)
{
    return SameValue<T>::check(a, b);
};

/// Assign the value of @p b into @p a using the type trait `AssignValue<T>` for type @p T.
template<typename T>
constexpr auto assignValue(T& a, const T& b)
{
    return AssignValue<T>::apply(a, b);
};

/// Return true if corresponding items in each tuple have the same value in the sense of `SameValue`.
template<typename Tuple1, typename Tuple2>
auto sameValues(const Tuple1& tuple1, const Tuple2& tuple2)
{
    using std::tuple_size_v;
    using std::get;
    constexpr auto N1 = tuple_size_v<Tuple1>;
    constexpr auto N2 = tuple_size_v<Tuple2>;
    static_assert(N1 == N2);
    bool res = true;
    For<N1>([&](auto i) constexpr {
        res = res && sameValue(get<i>(tuple1), get<i>(tuple2));
    });
    return res;
}

/// Assign the values of @p tuple2 into @p tuple1 in the sense of `AssignValue`.
template<typename Tuple1, typename Tuple2>
auto assignValues(Tuple1& tuple1, const Tuple2& tuple2)
{
    using std::tuple_size_v;
    using std::get;
    constexpr auto N1 = tuple_size_v<Tuple1>;
    constexpr auto N2 = tuple_size_v<Tuple2>;
    static_assert(N1 == N2);
    For<N1>([&](auto i) constexpr {
        assignValue(get<i>(tuple1), get<i>(tuple2));
    });
}

} // namespace detail

/// The class used to control memoization in the application.
class Memoization
{
public:
    /// Return true if memoization is currently enabled.
    static auto isEnabled() -> bool;

    /// Return true if memoization is currently disabled.
    static auto isDisabled() -> bool;

    /// Enable memoization optimization.
    static auto enable() -> void;

    /// Disable memoization optimization.
    static auto disable() -> void;

    /// Deleted default constructor.
    Memoization() = delete;
};

/// Return a memoized version of given function @p f.
template<typename Ret, typename... Args>
auto memoize(Fn<Ret(Args...)> f) -> Fn<Ret(Args...)>
{
    auto cache = std::make_shared<Map<Tuple<Args...>, Ret>>();
    return [=](Args... args) mutable -> Ret
    {
        if(Memoization::isDisabled())
            return f(args...);
        Tuple<Args...> t(args...);
        if(cache->find(t) == cache->end())
            (*cache)[t] = f(args...);
        return (*cache)[t];
    };
}

/// Return a memoized version of given function @p f.
template<typename Fun, EnableIf<!isFunction<Fun>>...>
auto memoize(Fun f)
{
    return memoize(asFunction(f));
}

/// Return a memoized version of given function @p f that caches only the arguments used in the last call.
template<typename Ret, typename... Args>
auto memoizeLast(Fn<Ret(Args...)> f) -> Fn<Ret(Args...)>
{
    Tuple<Decay<Args>...> cache;
    Ret result = Ret();
    auto firsttime = true;
    return [=](Args... args) mutable -> Ret
    {
        if(Memoization::isDisabled())
            return f(args...);
        if(detail::sameValues(cache, std::tie(args...)) && !firsttime)
            return Ret(result);
        detail::assignValues(cache, std::tie(args...));
        firsttime = false;
        return result = f(args...);
    };
}

/// Return a memoized version of given function @p f that caches only the arguments used in the last call.
template<typename Fun, EnableIf<!isFunction<Fun>>...>
auto memoizeLast(Fun f)
{
    return memoizeLast(asFunction(f));
}

/// Return a memoized version of given function @p f that caches only the arguments used in the last call.
template<typename Ret, typename RetRef, typename... Args>
auto memoizeLastUsingRef(Fn<void(RetRef, Args...)> f) -> Fn<void(RetRef, Args...)>
{
    Tuple<Decay<Args>...> cache;
    Ret result = Ret();
    auto firsttime = true;
    return [=](RetRef res, Args... args) mutable -> void
    {
        if(Memoization::isDisabled())
            f(res, args...);
        else if(detail::sameValues(cache, std::tie(args...)) && !firsttime)
            res = result;
        else
        {
            f(res, args...);
            result = res;
            detail::assignValues(cache, std::tie(args...));
        }
        firsttime = false;
    };
}

/// Return a memoized version of given function @p f that caches only the arguments used in the last call./// Return a memoized version of given function @p f that caches only the arguments used in the last call.
/// This overload is used when @p f is a lambda function or free function.
/// Use `memoizeLastUsingRef<Ret>(f)` to explicitly specify the `Ret` type.
template<typename Ret, typename Fun, EnableIf<!isFunction<Fun>>...>
auto memoizeLastUsingRef(Fun f)
{
    return memoizeLastUsingRef<Ret>(asFunction(f));
}

/// Return a memoized version of given function @p f that caches only the arguments used in the last call.
/// This overload assumes that `RetRef = Ret&`.
template<typename Ret, typename... Args>
auto memoizeLastUsingRef(Fn<void(Ret&, Args...)> f) -> Fn<void(Ret&, Args...)>
{
    return memoizeLastUsingRef<Ret, Ret&>(f);
}

/// Return a memoized version of given function @p f that caches only the arguments used in the last call.
/// This overload is used when @p f is a lambda function or free function.
/// Use `memoizeLastUsingRef(f)` to implicitly specify that `RetRef` is `Ret&`.
template<typename Fun, EnableIf<!isFunction<Fun>>...>
auto memoizeLastUsingRef(Fun f)
{
    return memoizeLastUsingRef(asFunction(f));
}

} // namespace Reaktoro
