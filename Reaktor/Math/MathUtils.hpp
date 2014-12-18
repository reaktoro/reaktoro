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

namespace Reaktor {

/// Determine the set of linearly independent columns in a matrix
/// @param A The matrix whose linearly independent columns must be found
/// @return The indices of the linearly independent columns
Indices linearlyIndependentCols(const Matrix& A);

/// Determine the set of linearly independent rows in a matrix
/// @param A The matrix whose linearly independent rows must be found
/// @return The indices of the linearly independent rows
Indices linearlyIndependentRows(const Matrix& A);

/// Determine the set of linearly independent columns in a matrix
/// @param[in] A The matrix whose linearly independent columns must be found
/// @param[out] B The matrix composed by linearly independent columns only
/// @return The indices of the linearly independent columns
Indices linearlyIndependentCols(const Matrix& A, Matrix& B);

/// Determine the set of linearly independent rows in a matrix
/// @param[in] A The matrix whose linearly independent rows must be found
/// @param[out] B The matrix composed by linearly independent rows only
/// @return The indices of the linearly independent rows
Indices linearlyIndependentRows(const Matrix& A, Matrix& B);

} // namespace Reaktor
