// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

/// A class that defines a Tridiagonal Matrix used on TransportSolver.
/// It stores data in a Eigen::VectorXd like, M = {a[0][0], a[0][1], a[0][2],
///                                                a[1][0], a[1][1], a[1][2],
///                                                a[2][0], a[2][1], a[2][2]}
class TridiagonalMatrix
{
public:
    TridiagonalMatrix() : TridiagonalMatrix(0) {}

    TridiagonalMatrix(Index size) : m_size(size), m_data(size * 3) {}

    auto size() const -> Index { return m_size; }

    auto data() -> VectorRef { return m_data; }

    auto data() const -> VectorConstRef { return m_data; }

    auto row(Index index) -> VectorRef { return m_data.segment(3 * index, 3); }

    auto row(Index index) const -> VectorConstRef { return m_data.segment(3 * index, 3); }

    auto a() -> VectorStridedRef { return Vector::Map(m_data.data() + 3, size() - 1, Eigen::InnerStride<3>()); }

    auto a() const -> VectorConstRef { return Vector::Map(m_data.data() + 3, size() - 1, Eigen::InnerStride<3>()); }

    auto b() -> VectorStridedRef { return Vector::Map(m_data.data() + 1, size(), Eigen::InnerStride<3>()); }

    auto b() const -> VectorConstRef { return Vector::Map(m_data.data() + 1, size(), Eigen::InnerStride<3>()); }

    auto c() -> VectorStridedRef { return Vector::Map(m_data.data() + 2, size() - 1, Eigen::InnerStride<3>()); }

    auto c() const -> VectorConstRef { return Vector::Map(m_data.data() + 2, size() - 1, Eigen::InnerStride<3>()); }

    auto resize(Index size) -> void;

    /// Factorize the tridiagonal matrix to help with linear system A x = b.
    auto factorize() -> void;

    /// Solve a linear system Ax = b with LU decomposition.
    auto solve(VectorRef x, VectorConstRef b) const -> void;

    /// Solve a linear system Ax = b with LU decomposition, using x as the unknown and it's.
    /// old values as the vector b.
    auto solve(VectorRef x) const -> void;

    operator Matrix() const;

private:
    /// The size of the tridiagonal matrix
    Index m_size;

    /// The coefficients
    Vector m_data;
};

} // namespace Reaktoro
