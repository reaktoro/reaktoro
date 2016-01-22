// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

#include "LU.hpp"

// Eigen includes
#include <Reaktoro/Eigen/Dense>

namespace Reaktoro {

LU::LU()
{}

LU::LU(const Matrix& A)
{
    compute(A);
}

LU::LU(const Matrix& A, const Vector& W)
{
    compute(A, W);
}

auto LU::empty() const -> bool
{
    return L.size();
}

auto LU::compute(const Matrix& A) -> void
{
    // The number of rows and cols of A
    const Index m = A.rows();
    const Index n = A.cols();
    const Index r = std::min(m, n);

    // Compute the full-pivoting LU of A
    Eigen::FullPivLU<Matrix> lu(A);

    // Set the rank of the formula matrix Ae
    rank = lu.rank();

    // Initialize the L, U, P, Q matrices so that P*A*Q = L*U
    L = lu.matrixLU().leftCols(r).triangularView<Eigen::UnitLower>();
    U = lu.matrixLU().topRows(r).triangularView<Eigen::Upper>();
    P = lu.permutationP();
    Q = lu.permutationQ();
}

auto LU::compute(const Matrix& A, const Vector& W) -> void
{
    // The number of rows and cols of A
    const Index m = A.rows();
    const Index n = A.cols();
    const Index r = std::min(m, n);

    // Initialize the weighted formula matrix
    const Matrix AW = A * diag(W);

    // Compute the full-pivoting LU of A
    Eigen::FullPivLU<Matrix> lu(tr(AW));

    // Set the rank of the matrix A
    rank = lu.rank();

    // Initialize the L, U, P, Q matrices so that P*A*Q = L*U
    L = tr(lu.matrixLU()).leftCols(r).triangularView<Eigen::Lower>();
    U = tr(lu.matrixLU()).triangularView<Eigen::UnitUpper>();
    P = lu.permutationQ().inverse();
    Q = lu.permutationP().inverse();

    // Correct the U matrix by unscaling it by weights
    U = U * Q.inverse() * diag(inv(W)) * Q;
}

auto LU::solve(const Vector& b) -> Vector
{
    const Index n = U.cols();
    Vector x = zeros(n);
    auto xx = x.segment(0, rank);
    xx = (P * b).segment(0, rank);
    xx = L.topLeftCorner(rank, rank).triangularView<Eigen::Lower>().solve(xx);
    xx = U.topLeftCorner(rank, rank).triangularView<Eigen::Upper>().solve(xx);
    x = Q * x;
    return x;
}

auto LU::trsolve(const Vector& b) -> Vector
{
    const Index m = L.rows();
    const auto& indices = Q.indices();
    Vector x(m);
    for(Index i = 0; i < m; ++i)
        x[i] = b[indices[i]];
    auto xx = x.segment(0, rank);
    xx = tr(U).topLeftCorner(rank, rank).triangularView<Eigen::Lower>().solve(xx);
    xx = tr(L).topLeftCorner(rank, rank).triangularView<Eigen::Upper>().solve(xx);
    x = P.inverse() * x;
    return x;
}

} // namespace Reaktoro
