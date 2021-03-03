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

template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
struct isArrayAux<Eigen::Array<Scalar, Rows, Cols, Options, MaxRows, MaxCols>> : std::true_type {};

template<typename Vec, int Size>
struct isArrayAux<Eigen::VectorBlock<Vec, Size>> : isArrayAux<Decay<Vec>> {};

template<typename Vec>
struct isArrayAux<Eigen::Ref<Vec>> : isArrayAux<Decay<Vec>> {};

/// Return the length of an argument (array size if array, 1 otherwise).
template<typename Arg>
constexpr auto length(const Arg& x)
{
    if constexpr(isArray<Arg>)
        return x.size();
    else return 1;
}

template<typename Arg, typename... Args>
constexpr auto length(const Arg& x, const Args&... xs)
{
    if constexpr(sizeof...(Args) > 0)
        return length(x) + length(xs...);
    else return length(x);
}

template<typename Array, typename Arg>
constexpr auto from(Array& array, std::size_t pos, const Arg& x)
{
    if constexpr(isArray<Arg>)
        array.segment(pos, x.size()) = x;
    else array[pos] = x;
}

template<typename Array, typename Arg, typename... Args>
constexpr auto from(Array& array, std::size_t pos, const Arg& x, const Args&... xs)
{
    from(array, pos, x);
    if constexpr(sizeof...(Args) > 0)
        from(array, pos + length(x), xs...);
}

template<typename Array, typename Arg>
constexpr auto to(const Array& array, std::size_t pos, Arg& x)
{
    if constexpr(isArray<Arg>)
        x = array.segment(pos, x.size());
    else x = array[pos];
}

template<typename Array, typename Arg, typename... Args>
constexpr auto to(const Array& array, std::size_t pos, Arg& x, Args&... xs)
{
    to(array, pos, x);
    if constexpr(sizeof...(Args) > 0)
        to(array, pos + length(x), xs...);
}

} // namespace detail

/// The class used to serialize/deserialize data using array.
template<typename T>
class ArrayStream
{
public:
    /// The array type used to store data.
    using ArrayType = ArrayX<Decay<T>>;

    /// Construct a default ArrayStream object.
    ArrayStream() = default;

    /// Construct an ArrayStream object with given list of scalars and/or arrays.
    template<typename Arg, typename... Args>
    ArrayStream(const Arg& x, const Args&... xs)
    : array(detail::length(x, xs...))
    {
        detail::from(array, 0, x, xs...);
    }

    /// Initialize this ArrayStream object with given list of scalars and/or arrays.
    template<typename Arg, typename... Args>
    auto from(const Arg& x, const Args&... xs)
    {
        array.resize(detail::length(x, xs...));
        detail::from(array, 0, x, xs...);
    }

    /// Transfer this ArrayStream object data to given list of scalars and/or arrays.
    template<typename Arg, typename... Args>
    auto to(Arg& x, Args&... xs) const
    {
        detail::to(array, 0, x, xs...);
    }

    auto data() const -> const ArrayType&
    {
        return array;
    }

    operator const ArrayType&() const
    {
        return array;
    }

private:
    /// The assembled array from a list of scalars and/or arrays.
    ArrayType array;
};

} // namespace Reaktoro
