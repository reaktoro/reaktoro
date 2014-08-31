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

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/Matrix.hpp>
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

/// Sets the specified entries of a vector with the entries of another
inline auto setRows(const Indices& irows, const Vector& values, Vector& vec) -> void
{
    for(unsigned i = 0; i < irows.size(); ++i)
        vec[irows[i]] = values[i];
}

/// Extracts the specified rows of the vector
inline auto rows(const Indices& irows, const Vector& vec) -> Vector
{
    Vector res(irows.size());

    for(unsigned i = 0; i < irows.size(); ++i)
        res[i] = vec[irows[i]];

    return res;
}

/// Extracts the specified rows of the matrix
inline auto rows(const Indices& irows, const Matrix& mat) -> Matrix
{
    Matrix res(irows.size(), mat.cols());

    for(unsigned i = 0; i < irows.size(); ++i)
        res.row(i) = mat.row(irows[i]);

    return res;
}

/// Extracts the specified columns of the matrix
inline auto cols(const Indices& icols, const Matrix& mat) -> Matrix
{
    Matrix res(mat.rows(), icols.size());

    for(unsigned j = 0; j < icols.size(); ++j)
        res.col(j) = mat.col(icols[j]);

    return res;
}

/// Extracts the specified rows and columns of the matrix
inline auto submatrix(const Indices& irows, const Indices& icols, const Matrix& mat) -> Matrix
{
    Matrix res(irows.size(), icols.size());

    for(unsigned i = 0; i < irows.size(); ++i)
        for(unsigned j = 0; j < icols.size(); ++j)
            res(i, j) = mat(irows[i], icols[j]);

    return res;
}

} // namespace Reaktor
