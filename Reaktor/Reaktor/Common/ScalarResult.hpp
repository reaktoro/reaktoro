// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

// Forward declarations
class VectorResultRow;

/**
 * Define a type that represents the result of a scalar function and its gradient and Hessian
 *
 * A @ref ScalarResult instance is a convenient way of expressing the result of a scalar function evaluation
 * that possibly includes its gradient and Hessian. See below for an example of its usage.
 *
 * @code
 * ScalarResult myScalarFunction(double x, double y)
 * {
 *     // Creates a ScalarResult instance with dimension 2
 *     ScalarResult res(0.0, zeros(2), zeros(2, 2));
 *
 *     // Sets the function result of the evaluation
 *     res.val = x * y;
 *
 *     // Sets the gradient result of the evaluation
 *     res.grad << y, x;
 *
 *     // Sets the Hessian result of the evaluation
 *     res.hessian << 0, 1,
 *                    1, 0;
 *
 *     return res;
 * }
 * @endcode
 *
 * @see VectorResult, Vector, Matrix
 */
class ScalarResult
{
public:
	/// Construct a default ScalarResult instance
	ScalarResult();

	/// Construct a ScalarResult instance
	/// @param dim The dimension of the domain of the scalar function
	ScalarResult(unsigned dim);

	/// Construct a ScalarResult instance
	/// @param val The scalar value of the evaluation of a scalar function
	/// @param grad The gradient of the evaluation of a scalar function
	ScalarResult(double val, const Vector& grad);

	/// Construct a ScalarResult instance from a row of a VectorResult instance
	/// @param res The VectorResultRow instance from which the ScalarResult instance is built
	ScalarResult(const VectorResultRow& res);

	/// The scalar value of the evaluation of a scalar function
	double func;

	/// The gradient of the evaluation of a scalar function
	Vector grad;
};

} // namespace Reaktor
