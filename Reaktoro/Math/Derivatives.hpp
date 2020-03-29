// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

/// Define a scalar function type
using ScalarFunction = std::function<double(VectorXdConstRef)>;

/// Define a vector function type
using VectorFunction = std::function<VectorXd(VectorXdConstRef)>;

/// Calculate the partial derivatives of a scalar function using a 1st-order forward finite difference scheme
///
/// @param f The scalar function as an instance of @ref ScalarFunction
/// @param x The point where the derivative is calculated
///
/// @return The partial derivatives of the scalar function
auto derivativeForward(const ScalarFunction& f, VectorXdConstRef x) -> VectorXd;

/// Calculate the partial derivatives of a scalar function using a 1st-order backward finite difference scheme
///
/// @param f The scalar function as an instance of @ref ScalarFunction
/// @param x The point where the derivative is calculated
///
/// @return The partial derivatives of the scalar function
auto derivativeBackward(const ScalarFunction& f, VectorXdConstRef x) -> VectorXd;

/// Calculate the partial derivatives of a scalar function using a 2nd-order central finite difference scheme
///
/// @param f The scalar function as an instance of @ref ScalarFunction
/// @param x The point where the derivative is calculated
///
/// @return The partial derivatives of the scalar function
auto derivativeCentral(const ScalarFunction& f, VectorXdConstRef x) -> VectorXd;

/// Calculate the partial derivatives of a vector function using a 1st-order forward finite difference scheme
///
/// @param f The scalar function as an instance of @ref ScalarFunction
/// @param x The point where the derivative is calculated
///
/// @return The partial derivatives of the scalar function
auto derivativeForward(const VectorFunction& f, VectorXdConstRef x) -> MatrixXd;

/// Calculate the partial derivatives of a vector function using a 1st-order backward finite difference scheme
///
/// @param f The scalar function as an instance of @ref ScalarFunction
/// @param x The point where the derivative is calculated
///
/// @return The partial derivatives of the scalar function
auto derivativeBackward(const VectorFunction& f, VectorXdConstRef x) -> MatrixXd;

/// Calculate the partial derivatives of a vector function using a 2nd-order central finite difference scheme
///
/// @param f The scalar function as an instance of @ref ScalarFunction
/// @param x The point where the derivative is calculated
///
/// @return The partial derivatives of the scalar function
auto derivativeCentral(const VectorFunction& f, VectorXdConstRef x) -> MatrixXd;

} // namespace Reaktoro
