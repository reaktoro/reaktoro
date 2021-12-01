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

#include "TridiagonalMatrix.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

namespace internal {

template<typename ReturnType, typename TridiagonalMatrixType>
auto row(TridiagonalMatrixType&& mat, Index index) -> ReturnType
{
    auto n = mat.size();
    auto data = mat.data();
    auto length = data.size();
    return (index == 0) ? data.segment(1, 2) : (index == n - 1) ?
        data.segment(length - 4, 2) : data.segment(3 * index, 3);
}

} // namespace internal

auto TridiagonalMatrix::resize(Index size) -> void
{
    m_size = size;
    m_data.conservativeResize(size * 3);
}

auto TridiagonalMatrix::factorize() -> void
{
    const Index n = size();

    auto prev = row(0).data(); // iterator to previous row
    auto curr = row(1).data(); // iterator to current row

    for(Index i = 1; i < n; ++i, prev = curr, curr += 3)
    {
        const auto& b_prev = prev[1]; // `b` value on the previous row
        const auto& c_prev = prev[2]; // `c` value on the previous row

        auto& a_curr = curr[0]; // `a` value on the current row
        auto& b_curr = curr[1]; // `b` value on the current row

        a_curr /= b_prev; // update the a-diagonal in the tridiagonal matrix
        b_curr -= a_curr * c_prev; // update the b-diagonal in the tridiagonal matrix
    }
}

auto TridiagonalMatrix::solve(ArrayXrRef x, ArrayXrConstRef d) const -> void
{
    const Index n = size();

    auto curr = row(1).data(); // iterator to current row

    //-------------------------------------------------------------------------
    // Perform the forward solve with the L factor of the LU factorization
    //-------------------------------------------------------------------------
    x[0] = d[0];

    for(Index i = 1; i < n; ++i, curr += 3)
    {
        const auto& a = curr[0]; // `a` value on the current row

        x[i] = d[i] - a * x[i - 1];
    }

    curr -= 3; // step back so that curr points to the last row
    const auto& bn = curr[1]; // `b` value on the last row
    curr -= 3; // step back so that curr points to the second to last row

    //-------------------------------------------------------------------------
    // Perform the backward solve with the U factor of the LU factorization
    //-------------------------------------------------------------------------
    x[n - 1] /= bn;

    for(Index i = 2; i <= n; ++i, curr -= 3)
    {
        const auto& k = n - i; // the index of the current row
        const auto& b = curr[1]; // `b` value on the current row
        const auto& c = curr[2]; // `c` value on the current row

        x[k] = (x[k] - c * x[k + 1])/b;
    }
}

auto TridiagonalMatrix::solve(ArrayXrRef x) const -> void
{
    solve(x, x);
}

} // namespace Reaktoro
