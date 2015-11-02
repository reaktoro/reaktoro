// Reaktoro is a C++ library for computational reaction modelling.
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

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

/// Determine the set of linearly independent columns in a matrix using a column pivoting QR algorithm.
/// @param A The matrix whose linearly independent columns must be found
/// @return The indices of the linearly independent columns
auto linearlyIndependentCols(const Matrix& A) -> Indices;

/// Determine the set of linearly independent columns in a matrix using partial pivoting LU algorithm.
/// @param A The matrix whose linearly independent columns must be found
/// @return The indices of the linearly independent columns
auto linearlyIndependentColsPartialPivLU(const Matrix& A) -> Indices;

/// Determine the set of linearly independent rows in a matrix.
/// @param A The matrix whose linearly independent rows must be found
/// @return The indices of the linearly independent rows
auto linearlyIndependentRows(const Matrix& A) -> Indices;

/// Determine the set of linearly independent columns in a matrix.
/// @param[in] A The matrix whose linearly independent columns must be found
/// @param[out] B The matrix composed by linearly independent columns only
/// @return The indices of the linearly independent columns
auto linearlyIndependentCols(const Matrix& A, Matrix& B) -> Indices;

/// Determine the set of linearly independent rows in a matrix.
/// @param[in] A The matrix whose linearly independent rows must be found
/// @param[out] B The matrix composed by linearly independent rows only
/// @return The indices of the linearly independent rows
auto linearlyIndependentRows(const Matrix& A, Matrix& B) -> Indices;

/// Calculate the inverse of `A + D` where `inv(A)` is already known and `D` is a diagonal matrix.
/// @param invA[in,out] The inverse of the matrix `A` and the final inverse of `A + D`
/// @param D The diagonal matrix `D`
auto inverseShermanMorrison(const Matrix& invA, const Vector& D) -> Matrix;

/// Calculate the LU decomposition of a square or non-square matrix using partial pivoting.
/// @param A The matrix to be decomposed.
/// @param L The lower triangular matrix `L` of the LU decomposition `PAQ = LU`.
/// @param U The upper triangular matrix `U` of the LU decomposition `PAQ = LU`.
/// @param P The permutation matrix `P` of the LU decomposition `PAQ = LU`.
/// @param Q The permutation matrix `Q` of the LU decomposition `PAQ = LU`.
auto lu(const Matrix& A, Matrix& L, Matrix& U, PermutationMatrix& P, PermutationMatrix& Q) -> void;

} // namespace Reaktoro
