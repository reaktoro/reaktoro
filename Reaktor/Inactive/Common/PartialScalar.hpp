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
#include <functional>

// Reaktor includes
#include <Reaktor/Common/PartialUtils.hpp>
#include <Reaktor/Common/Matrix.hpp>
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

/**
 * Defines a type that represents the result of a scalar function and its gradient and Hessian
 *
 * A @ref PartialScalar instance is a convenient way of expressing the result of a scalar function evaluation
 * that possibly includes its gradient and Hessian. See below for an example of its usage.
 *
 * @code
 * PartialScalar myScalarFunction(double x, double y)
 * {
 *     // Creates a PartialScalar instance with dimension 2
 *     PartialScalar res = partialScalar(0.0, zeros(2), zeros(2, 2));
 *
 *     // Sets the function result of the evaluation
 *     func(res) = x * y;
 *
 *     // Sets the gradient result of the evaluation
 *     grad(res) << y, x;
 *
 *     // Sets the Hessian result of the evaluation
 *     hessian(res) << 0, 1,
 *                     1, 0;
 *
 *     return res;
 * }
 * @endcode
 *
 * In the code above, a scalar function @c myScalarFunction is defined. The utility function
 * @ref partialScalar is used to create a @ref PartialScalar instance, while the functions @ref func,
 * @ref grad and @ref hessian are used to access the function, gradient and Hessian results of the
 * function evaluation.
 *
 * @see partialScalar, func, grad, hessian
 */
using PartialScalar = std::tuple<double, Vector, Matrix>;

/**
 * Creates a @ref PartialScalar instance
 * @param val The value of the function evaluation
 * @return A @ref PartialScalar instance with uninitialized gradient and Hessian
 */
template<typename Value>
inline auto partialScalar(Value&& val) -> PartialScalar
{
    return PartialScalar{std::forward<Value>(val), Vector(), Matrix()};
}

/**
 * Creates a @ref PartialScalar instance
 * @param val The value of the function evaluation
 * @param grad The gradient of the function evaluation
 * @return A @ref PartialScalar instance with uninitialized Hessian
 */
template<typename Value, typename Grad>
inline auto partialScalar(Value&& val, Grad&& grad) -> PartialScalar
{
    return PartialScalar{std::forward<Value>(val), std::forward<Grad>(grad), Matrix()};
}

/**
 * Creates a @ref PartialScalar instance
 * @param val The value of the function evaluation
 * @param grad The gradient of the function evaluation
 * @param hessian The Hessian of the function evaluation
 * @return A @ref PartialScalar instance
 */
template<typename Value, typename Grad, typename Hessian>
inline auto partialScalar(Value&& val, Grad&& grad, Hessian&& hessian) -> PartialScalar
{
    return PartialScalar{std::forward<Value>(val), std::forward<Grad>(grad), std::forward<Hessian>(hessian)};
}

} // namespace Reaktor
