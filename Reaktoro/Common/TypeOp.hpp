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

// C++ includes
#include <type_traits>

namespace Reaktoro {

namespace detail {

template<typename T>
struct TypeOpRef
{
    using type = T&;
};

template<typename T>
struct TypeOpConstRef
{
    using type = T const&;
};

template<typename T>
struct TypeOpIdentity
{
    using type = T;
};

} // namespace detail

template<typename T>
using TypeOpRef = typename detail::TypeOpRef<std::decay_t<T>>::type;

template<typename T>
using TypeOpConstRef = typename detail::TypeOpConstRef<std::decay_t<T>>::type;

template<typename T>
using TypeOpIdentity = typename detail::TypeOpIdentity<std::decay_t<T>>::type;

} // namespace Reaktoro
