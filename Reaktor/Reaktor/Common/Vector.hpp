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

/// Define an alias to the vector type of the Armadillo library
using Vector = arma::vec;

/// Define a type for a contiguous view of a vector
using VectorView = arma::subview_col<double>;

/// Define a type for a view of a single row of a vector
using VectorRow = arma::subview_col<double>;

/// Define a type for a non-contiguous view of a vector
using SubVector = arma::subview_elem1<double, arma::Mat<unsigned>>;

/// Return an expression of a zero vector
/// @param rows The number of rows
/// @return The expression of a zero vector
inline auto zeros(unsigned rows) -> decltype(arma::zeros<arma::vec>(rows))
{
	return arma::zeros<arma::vec>(rows);
}

/// Return an expression of a vector with entries equal to one
/// @param rows The number of rows
/// @return The expression of a vector with entries equal to one
inline auto ones(unsigned rows) -> decltype(arma::ones<arma::vec>(rows))
{
    return arma::ones<arma::vec>(rows);
}

} /* namespace Reaktor */
