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
#include <Reaktoro/Common/ArraySerialization.hpp>

namespace Reaktoro {

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
    template<typename... Args>
    ArrayStream(const Args&... args)
    {
        ArraySerialization::resize(array, args...);
        ArraySerialization::serialize(array, args...);
    }

    /// Initialize this ArrayStream object with given list of scalars and/or arrays.
    template<typename... Args>
    auto from(const Args&... args)
    {
        ArraySerialization::resize(array, args...);
        ArraySerialization::serialize(array, args...);
    }

    /// Transfer this ArrayStream object data to given list of scalars and/or arrays.
    template<typename... Args>
    auto to(Args&... args) const
    {
        ArraySerialization::deserialize(array, args...);
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
