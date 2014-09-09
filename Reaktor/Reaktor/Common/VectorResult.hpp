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
#include <Reaktor/Common/Matrix.hpp>
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

// Forward declarations
class ScalarResult;
class VectorResultRow;

/// Define a type that represents the result of a vector-valued function and its gradient and Hessian
///
/// A VectorResult instance is a convenient way of expressing the result of a vector-valued function
/// evaluation that possibly includes its gradient and Hessian. See below for an example of its usage.
///
/// @code
/// VectorResult myVectorFunction(double x, double y)
/// {
///     // Creates a VectorResult instance with dimension 2
///     VectorResult res(zeros(2), zeros(2, 2));
///
///     // Sets the function result of the evaluation
///     res.func << x + y, x - y;
///
///     // Sets the gradient result of the evaluation
///     res.grad << 1,  1,
///                 1, -1;
///
///     return res;
/// }
/// @endcode
///
/// @see ScalarResult
class VectorResult
{
public:
	/// Construct a default VectorResult instance
	VectorResult();

	/// Construct a VectorResult instance
	/// @param dim The dimension of the domain and range of the vector function
	VectorResult(unsigned dim);

	/// Construct a VectorResult instance
	/// @param domain The dimension of the domain of the vector function
	/// @param range The dimension of the range of the vector function
	VectorResult(unsigned domain, unsigned range);

	/// Construct a VectorResult instance
	/// @param val The vector value of the evaluation of a vector function
	/// @param grad The gradient of the evaluation of a vector function
	VectorResult(const Vector& val, const Matrix& grad);

	auto row(const Index& i) -> VectorResultRow;

	auto row(const Index& i) const -> const VectorResultRow;

	/// The vector value of the evaluation of a vector function
	Vector func;

	/// The gradient of the evaluation of a vector function
	Matrix grad;
};

class VectorResultRow
{
public:
	/// Construct a default VectorResultRow instance
	VectorResultRow();

	/// Construct a VectorResultRow instance
	/// @param val The vector value of the evaluation of a vector function
	/// @param grad The gradient of the evaluation of a vector function
	VectorResultRow(const VectorRow& val_row, const MatrixRow& grad_row);

	auto operator=(const ScalarResult& res) -> VectorResultRow&;

	/// The reference to the row of the vector value of the evaluation of a vector function
	VectorRow func;

	/// The reference to the row of the gradient of the evaluation of a vector function
	MatrixRow grad;
};

} // namespace Reaktor
