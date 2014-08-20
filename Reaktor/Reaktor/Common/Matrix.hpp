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

// Armadillo includes
#include <armadillo>

namespace Reaktor {

/// Define an alias to the matrix type of the Armadillo library
using Matrix = arma::mat;

/// Define an alias to the submatrix type of the Armadillo library
using SubMatrix = arma::subview<double>;

/// Define an alias to the matrix-row type of the Armadillo library
using MatrixRow = arma::subview_row<double>;

/// Define an alias to the matrix-col type of the Armadillo library
using MatrixCol = arma::subview_col<double>;

/// Return an expression of a zero matrix
/// @param rows The number of rows
/// @param cols The number of columns
/// @return The expression of a zero matrix
inline auto zeros(unsigned rows, unsigned cols) -> decltype(arma::zeros<arma::mat>(rows, cols))
{
    return arma::zeros<arma::mat>(rows, cols);
}

/// Return an expression of a matrix with entries equal to one
/// @param rows The number of rows
/// @param cols The number of columns
/// @return The expression of a matrix with entries equal to one
inline auto ones(unsigned rows, unsigned cols) -> decltype(arma::ones<arma::mat>(rows, cols))
{
    return arma::ones<arma::mat>(rows, cols);
}

} /* namespace Reaktor */
