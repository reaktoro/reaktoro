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
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

inline auto convert(const Indices& indices) -> arma::uvec
{
    arma::uvec vec(indices.size());
    std::copy(indices.begin(), indices.end(), vec.memptr());
    return vec;
}

/// Set the specified entries of a vector with the entries of another
inline auto setRows(const Indices& irows, const Vector& values, Vector& vec) -> void
{
    for(unsigned i = 0; i < irows.size(); ++i)
        vec[irows[i]] = values[i];
}

inline auto rows(Vector& v, unsigned start, unsigned size) -> decltype(v.rows(start, start+size-1))
{
    return v.rows(start, start+size-1);
}

inline auto rows(const Vector& v, unsigned start, unsigned size) -> Vector
{
    return v.rows(start, start+size-1);
}

inline auto rows(Matrix& m, const Indices& indices) -> decltype(m.rows(convert(indices)))
{
    return m.rows(convert(indices));
}

inline auto cols(Matrix& m, const Indices& indices) -> decltype(m.cols(convert(indices)))
{
    return m.cols(convert(indices));
}

inline auto rows(const Matrix& m, const Indices& indices) -> Matrix
{
    return m.rows(convert(indices));
}

inline auto cols(const Matrix& m, const Indices& indices) -> Matrix
{
    return m.cols(convert(indices));
}

/// Extract the specified rows and columns of the matrix
inline auto submatrix(const Matrix& m, const Indices& irows, const Indices& icols) -> Matrix
{
    return m.submat(convert(irows), convert(icols));
}

} // namespace Reaktor
