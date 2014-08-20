/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

// Eigen includes
#include <Eigen/Core>

namespace Reaktor {

/**
 * Defines an alias to the matrix type of the Eigen library
 */
using Matrix = Eigen::MatrixXd;

/**
 * Returns an expression of a zero matrix
 * @param rows The number of rows
 * @param cols The number of columns
 * @return The expression of a zero matrix
 */
inline auto zeros(unsigned rows, unsigned cols) -> decltype(Matrix::Zero(rows, cols))
{
    return Matrix::Zero(rows, cols);
}

/**
 * Returns an expression of a matrix with entries equal to one
 * @param rows The number of rows
 * @param cols The number of columns
 * @return The expression of the matrix with entries equal to one
 */
inline auto ones(unsigned rows, unsigned cols) -> decltype(Matrix::Ones(rows, cols))
{
    return Matrix::Ones(rows, cols);
}

} /* namespace Reaktor */
