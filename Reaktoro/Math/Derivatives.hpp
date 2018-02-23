// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

/// Define a scalar function type
using ScalarFunction = std::function<double(VectorConstRef)>;

/// Define a vector function type
using VectorFunction = std::function<Vector(VectorConstRef)>;

/// Calculate the partial derivatives of a scalar function using a 1st-order forward finite difference scheme
///
/// @param f The scalar function as an instance of @ref ScalarFunction
/// @param x The point where the derivative is calculated
///
/// @return The partial derivatives of the scalar function
auto derivativeForward(const ScalarFunction& f, VectorConstRef x) -> Vector;

/// Calculate the partial derivatives of a scalar function using a 1st-order backward finite difference scheme
///
/// @param f The scalar function as an instance of @ref ScalarFunction
/// @param x The point where the derivative is calculated
///
/// @return The partial derivatives of the scalar function
auto derivativeBackward(const ScalarFunction& f, VectorConstRef x) -> Vector;

/// Calculate the partial derivatives of a scalar function using a 2nd-order central finite difference scheme
///
/// @param f The scalar function as an instance of @ref ScalarFunction
/// @param x The point where the derivative is calculated
///
/// @return The partial derivatives of the scalar function
auto derivativeCentral(const ScalarFunction& f, VectorConstRef x) -> Vector;

/// Calculate the partial derivatives of a vector function using a 1st-order forward finite difference scheme
///
/// @param f The scalar function as an instance of @ref ScalarFunction
/// @param x The point where the derivative is calculated
///
/// @return The partial derivatives of the scalar function
auto derivativeForward(const VectorFunction& f, VectorConstRef x) -> Matrix;

/// Calculate the partial derivatives of a vector function using a 1st-order backward finite difference scheme
///
/// @param f The scalar function as an instance of @ref ScalarFunction
/// @param x The point where the derivative is calculated
///
/// @return The partial derivatives of the scalar function
auto derivativeBackward(const VectorFunction& f, VectorConstRef x) -> Matrix;

/// Calculate the partial derivatives of a vector function using a 2nd-order central finite difference scheme
///
/// @param f The scalar function as an instance of @ref ScalarFunction
/// @param x The point where the derivative is calculated
///
/// @return The partial derivatives of the scalar function
auto derivativeCentral(const VectorFunction& f, VectorConstRef x) -> Matrix;

} // namespace Reaktoro
