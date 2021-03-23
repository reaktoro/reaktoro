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

// Reaktoro includes
#include <Reaktoro/Common/TraitsUtils.hpp>
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {
namespace detail {

template<typename T>
struct isArrayAux;

template<typename T>
constexpr auto isArray = isArrayAux<Decay<T>>::value;

template<typename T>
struct isArrayAux : std::false_type {};

template<typename Scalar, int Rows, int Options, int MaxRows>
struct isArrayAux<Eigen::Array<Scalar, Rows, 1, Options, MaxRows, 1>> : std::true_type {};

template<typename Scalar, int Cols, int Options, int MaxCols>
struct isArrayAux<Eigen::Array<Scalar, 1, Cols, Options, 1, MaxCols>> : std::true_type {};

template<typename Scalar, int Rows, int Options, int MaxRows>
struct isArrayAux<Eigen::Matrix<Scalar, Rows, 1, Options, MaxRows, 1>> : std::true_type {};

template<typename Scalar, int Cols, int Options, int MaxCols>
struct isArrayAux<Eigen::Matrix<Scalar, 1, Cols, Options, 1, MaxCols>> : std::true_type {};

template<typename Vec, int Size>
struct isArrayAux<Eigen::VectorBlock<Vec, Size>> : isArrayAux<Decay<Vec>> {};

template<typename Vec>
struct isArrayAux<Eigen::Ref<Vec>> : isArrayAux<Decay<Vec>> {};

/// Return the length of an argument (array size if array, 1 otherwise).
template<typename Arg>
constexpr auto length(const Arg& x) -> std::size_t
{
    if constexpr(isArray<Arg>)
        return x.size();
    else return 1;
}

template<typename Arg, typename... Args>
constexpr auto length(const Arg& x, const Args&... xs) -> std::size_t
{
    if constexpr(sizeof...(Args) > 0)
        return length(x) + length(xs...);
    else return length(x);
}

template<typename Array, typename Arg>
constexpr auto serialize(Array& array, std::size_t pos, const Arg& x) -> void
{
    using U = Decay<decltype(array[0])>;

    if constexpr(isArray<Arg>)
    {
        assert(array.size() >= pos + x.size());
        array.segment(pos, x.size()) = x.template cast<U>();
    }
    else
    {
        assert(array.size() >= pos + 1);
        array[pos] = static_cast<U>(x);
    }
}

template<typename Array, typename Arg, typename... Args>
constexpr auto serialize(Array& array, std::size_t pos, const Arg& x, const Args&... xs) -> void
{
    serialize(array, pos, x);
    if constexpr(sizeof...(Args) > 0)
        serialize(array, pos + length(x), xs...);
}

template<typename Array, typename Arg>
constexpr auto deserialize(const Array& array, std::size_t pos, Arg& x) -> void
{
    if constexpr(isArray<Arg>)
    {
        assert(array.size() >= pos + x.size());
        x = array.segment(pos, x.size());
    }
    else
    {
        assert(array.size() >= pos + 1);
        x = array[pos];
    }
}

template<typename Array, typename Arg, typename... Args>
constexpr auto deserialize(const Array& array, std::size_t pos, Arg& x, Args&... xs) -> void
{
    deserialize(array, pos, x);
    if constexpr(sizeof...(Args) > 0)
        deserialize(array, pos + length(x), xs...);
}

} // namespace detail

/// The class implementing methods to serialize/deserialize data into/from arrays.
struct ArraySerialization
{
    /// Resize @p array so that it can store the data in @p args.
    template<typename Array, typename... Args>
    constexpr static auto resize(Array& array, const Args&... args) -> void
    {
        const auto size = detail::length(args...);
        array.resize(size);
    }

    /// Copy the data in @p args into @p array.
    template<typename Array, typename... Args>
    constexpr static auto serialize(Array& array, const Args&... args) -> void
    {
        detail::serialize(array, 0, args...);
    }

    /// Copy the data in @p array to @p args.
    template<typename Array, typename... Args>
    constexpr static auto deserialize(const Array& array, Args&... args) -> void
    {
        detail::deserialize(array, 0, args...);
    }
};

} // namespace Reaktoro
