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
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

/// Determine the set of linearly independent columns in a matrix using a column pivoting QR algorithm.
/// @param A The matrix whose linearly independent columns must be found
/// @return The indices of the linearly independent columns
auto linearlyIndependentCols(MatrixXdConstRef A) -> Indices;

/// Determine the set of linearly independent rows in a matrix.
/// @param A The matrix whose linearly independent rows must be found
/// @return The indices of the linearly independent rows
auto linearlyIndependentRows(MatrixXdConstRef A) -> Indices;

/// Determine the set of linearly independent columns in a matrix.
/// @param[in] A The matrix whose linearly independent columns must be found
/// @param[out] B The matrix composed by linearly independent columns only
/// @return The indices of the linearly independent columns
auto linearlyIndependentCols(MatrixXdConstRef A, MatrixXdRef B) -> Indices;

/// Determine the set of linearly independent rows in a matrix.
/// @param[in] A The matrix whose linearly independent rows must be found
/// @param[out] B The matrix composed by linearly independent rows only
/// @return The indices of the linearly independent rows
auto linearlyIndependentRows(MatrixXdConstRef A, MatrixXdRef B) -> Indices;

/// Calculate the inverse of `A + D` where `inv(A)` is already known and `D` is a diagonal matrix.
/// @param invA[in,out] The inverse of the matrix `A` and the final inverse of `A + D`
/// @param D The diagonal matrix `D`
auto inverseShermanMorrison(MatrixXdConstRef invA, VectorXdConstRef D) -> MatrixXd;

/// Calculates the rational number that approximates a given real number.
/// The algorithm is based on Farey sequence as shown
/// [here](http://www.johndcook.com/blog/2010/10/20/best-rational-approximation/).
/// @param x The real number.
/// @param maxden The maximum denominator.
/// @return A tuple containing the numerator and denominator.
auto rationalize(double x, unsigned maxden) -> std::tuple<long, long>;

/// Clean an array that is known to have rational numbers from round-off errors.
/// @param vals[in,out] The values to be cleaned
/// @param maxden The maximum known denominator in the array with rational numbers
auto cleanRationalNumbers(double* vals, long size, long maxden = 6) -> void;

/// Clean a matrix that is known to have rational numbers from round-off errors.
/// @param A[in,out] The matrix to be cleaned
/// @param maxden The maximum known denominator in the matrix with rational numbers
auto cleanRationalNumbers(MatrixXdRef A, long maxden = 6) -> void;

/// Return the dot product `s + dot(x, y)` of two vectors with triple-precision.
auto dot3p(VectorXdConstRef x, VectorXdConstRef y, double s) -> double;

/// Return the residual of the equation `A*x - b` with triple-precision.
auto residual3p(MatrixXdConstRef A, VectorXdConstRef x, VectorXdConstRef b) -> VectorXd;

} // namespace Reaktoro
