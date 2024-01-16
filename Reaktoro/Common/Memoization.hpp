// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

/// Used to enable function arguments of a type `T` to be cached in a memoized version of the function.
template<typename T>
struct MemoizationTraits
{
    /// The type `T` without const, reference, and pointer specifiers.
    using Type = Decay<T>;

    /// The type of the object constructed with a given object of type `T`.
    /// Consider three cases for type `T` and the respective default type for `CacheType` using `Type = Decay<T>`:
    ///
    /// **Case 1:** `T` is `const real&` and `CacheType` is `real`.
    /// **Case 2:** `T` is `ArrayXrConstRef` and `CacheType` is also `ArrayXrConstRef`.
    ///
    /// In Case 1, the default type for `CacheType` is suitable, since `real`
    /// implements a default constructor, which will be used to construct a
    /// cache object for arguments of this type.
    ///
    /// In Case 2, `CacheType` is not a suitable type for a cache object,
    /// because `ArrayXrConstRef` does not implement a default constructor.
    /// In this case, `CacheType` should be defined to be `ArrayXr` instead.
    using CacheType = Type;

    /// Check if `a`, of type `CacheType`, is equal to `b`, of type `T`.
    /// Redefine this function in a template specialization for type `T` if the
    /// default implementation using `operator==(const CacheType&, const T&)` is not adequate.
    /// For example, `operator==` for `Eigen::ArrayBase<Derived>` objects
    /// return an expression that needs to be evaluated with methods `all()` or
    /// `any()`. Thus, in order for memoization of functions with this array
    /// type to compile, a template specialization is needed for `Eigen::ArrayBase<Derived>>`.
    static auto equal(const CacheType& a, const Type& b)
    {
        return a == b;
    }

    /// Assign the value of `b`, of type `T`, into `a`, of type `CacheType`.
    /// Redefine this function in a template specialization for type `T` if the
    /// default implementation using `CacheType::operator=` is not adequate.
    static auto assign(CacheType& a, const Type& b)
    {
        a = b;
    }
};

namespace detail {

/// Return true if `a` and `b` have the same value using `MemoizationTraits::equal`.
template<typename CacheType, typename T>
constexpr auto sameValue(const CacheType& a, const T& b)
{
    return MemoizationTraits<Decay<T>>::equal(a, b);
};

/// Assign the value of `b` into `a` using `MemoizationTraits::assign`.
template<typename CacheType, typename T>
constexpr auto assignValue(CacheType& a, const T& b)
{
    return MemoizationTraits<Decay<T>>::assign(a, b);
};

/// Return true if corresponding items in each tuple have the same value.
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

/// Assign the values of `tuple2` into `tuple1`.
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

/// Used to get the type used as cache type in memoization of functions.
template<typename T>
using CacheType = typename MemoizationTraits<Decay<T>>::CacheType;

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

/// Return a memoized version of given function `f`.
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

/// Return a memoized version of given function `f`.
template<typename Fun, Requires<!isFunction<Fun>> = true>
auto memoize(Fun f)
{
    return memoize(asFunction(f));
}

/// Return a memoized version of given function `f` that caches only the arguments used in the last call.
template<typename Ret, typename... Args>
auto memoizeLast(Fn<Ret(Args...)> f) -> Fn<Ret(Args...)>
{
    Tuple<detail::CacheType<Args>...> cache;
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

/// Return a memoized version of given function `f` that caches only the arguments used in the last call.
template<typename Fun, Requires<!isFunction<Fun>> = true>
auto memoizeLast(Fun f)
{
    return memoizeLast(asFunction(f));
}

/// Return a memoized version of given function `f` that caches only the arguments used in the last call.
template<typename Ret, typename RetRef, typename... Args>
auto memoizeLastUsingRef(Fn<void(RetRef, Args...)> f) -> Fn<void(RetRef, Args...)>
{
    Tuple<detail::CacheType<Args>...> cache;
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

/// Return a memoized version of given function `f` that caches only the arguments used in the last call./// Return a memoized version of given function `f` that caches only the arguments used in the last call.
/// This overload is used when `f` is a lambda function or free function.
/// Use `memoizeLastUsingRef<Ret>(f)` to explicitly specify the `Ret` type.
template<typename Ret, typename Fun, Requires<!isFunction<Fun>> = true>
auto memoizeLastUsingRef(Fun f)
{
    return memoizeLastUsingRef<Ret>(asFunction(f));
}

/// Return a memoized version of given function `f` that caches only the arguments used in the last call.
/// This overload assumes that `RetRef = Ret&`.
template<typename Ret, typename... Args>
auto memoizeLastUsingRef(Fn<void(Ret&, Args...)> f) -> Fn<void(Ret&, Args...)>
{
    return memoizeLastUsingRef<Ret, Ret&>(f);
}

/// Return a memoized version of given function `f` that caches only the arguments used in the last call.
/// This overload is used when `f` is a lambda function or free function.
/// Use `memoizeLastUsingRef(f)` to implicitly specify that `RetRef` is `Ret&`.
template<typename Fun, Requires<!isFunction<Fun>> = true>
auto memoizeLastUsingRef(Fun f)
{
    return memoizeLastUsingRef(asFunction(f));
}

} // namespace Reaktoro
