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
using Cube = arma::cube;

/// Define an alias to the subcube type of the Armadillo library
using SubCube = arma::subview_cube<double>;

/// Return an expression of a zero cube
/// @param rows The number of rows
/// @param cols The number of columns
/// @param slices The number of slices
/// @return The expression of a zero cube
inline auto zeros(unsigned rows, unsigned cols, unsigned slices) -> decltype(arma::zeros<arma::cube>(rows, cols, slices))
{
    return arma::zeros<arma::cube>(rows, cols, slices);
}

/// Return an expression of a vector with entries equal to one
/// @param rows The number of rows
/// @param cols The number of columns
/// @param slices The number of slices
/// @return The expression of a cube with entries equal to one
inline auto ones(unsigned rows, unsigned cols, unsigned slices) -> decltype(arma::ones<arma::cube>(rows, cols, slices))
{
    return arma::ones<arma::cube>(rows, cols, slices);
}

} /* namespace Reaktor */


