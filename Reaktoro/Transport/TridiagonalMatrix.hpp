// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

/// A class that defines a Tridiagonal Matrix used on TransportSolver.
/// It stores data in a Eigen::ArrayXrXd like, M = {a[0][0], a[0][1], a[0][2],
///                                                a[1][0], a[1][1], a[1][2],
///                                                a[2][0], a[2][1], a[2][2]}
class TridiagonalMatrix
{
public:
    /// Construct a default TridiagonalMatrix instance.
    TridiagonalMatrix() : TridiagonalMatrix(0) {}

    /// Construct a TridiagonalMatrix instance of a given size.
    TridiagonalMatrix(Index size) : m_size(size), m_data(size * 3) {}

    /// Return the size of the matrix.
    auto size() const -> Index { return m_size; }

    /// Return the coefficients of the TridiagonalMatrix instance
    auto data() -> ArrayXrRef { return m_data; }

    /// Return the const coefficients of the TridiagonalMatrix instance
    auto data() const -> ArrayXrConstRef { return m_data; }

    /// Return the row with a given index
    auto row(Index index) -> ArrayXrRef { return m_data.segment(3 * index, 3); }

    /// Return the const row with a given index
    auto row(Index index) const -> ArrayXrConstRef { return m_data.segment(3 * index, 3); }

    /// Return an entrance that corresponds to the a-block
    auto a() -> ArrayXrStridedRef { return ArrayXr::Map(m_data.data() + 3, size() - 1, Eigen::InnerStride<3>()); }

    /// Return a constant entrance that corresponds to the a-block
    auto a() const -> ArrayXrConstRef { return ArrayXr::Map(m_data.data() + 3, size() - 1, Eigen::InnerStride<3>()); }

    /// Return an entrance that corresponds to the b-block
    auto b() -> ArrayXrStridedRef { return ArrayXr::Map(m_data.data() + 1, size(), Eigen::InnerStride<3>()); }

    /// Return a constant entrance that corresponds to the a-block
    auto b() const -> ArrayXrConstRef { return ArrayXr::Map(m_data.data() + 1, size(), Eigen::InnerStride<3>()); }

    /// Return an entrance that corresponds to the c-block
    auto c() -> ArrayXrStridedRef { return ArrayXr::Map(m_data.data() + 2, size() - 1, Eigen::InnerStride<3>()); }

    /// Return a const entrance that corresponds to the c-block
    auto c() const -> ArrayXrConstRef { return ArrayXr::Map(m_data.data() + 2, size() - 1, Eigen::InnerStride<3>()); }

    /// Resize function
    auto resize(Index size) -> void;

    /// Factorize the tridiagonal matrix to help with linear system A x = b.
    auto factorize() -> void;

    /// Solve a linear system Ax = b with LU decomposition.
    auto solve(ArrayXrRef x, ArrayXrConstRef b) const -> void;

    /// Solve a linear system Ax = b with LU decomposition, using x as the unknown and it's.
    /// old values as the vector b.
    auto solve(ArrayXrRef x) const -> void;

private:
    /// The size of the tridiagonal matrix.
    Index m_size;

    /// The coefficients.
    ArrayXr m_data;
};

} // namespace Reaktoro
