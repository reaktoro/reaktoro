/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

// C++ includes
#include <tuple>

namespace Reaktor {

/**
 * Access the function component of a @ref ScalarResult or @ref VectorResult instance
 */
template<typename TupleType>
inline auto func(TupleType&& val) -> decltype(std::get<0>(val))
{
    return std::get<0>(val);
}

/**
 * Access the gradient component of a @ref ScalarResult or @ref VectorResult instance
 */
template<typename TupleType>
inline auto grad(TupleType&& val) -> decltype(std::get<1>(val))
{
    return std::get<1>(val);
}

/**
 * Access the Hessian component of a @ref ScalarResult or @ref VectorResult instance
 */
template<typename TupleType>
inline auto hessian(TupleType&& val) -> decltype(std::get<2>(val))
{
    return std::get<2>(val);
}


} // namespace Reaktor
