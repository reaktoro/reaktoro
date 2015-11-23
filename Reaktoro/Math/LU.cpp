// Reaktoro is a C++ library for computational reaction modelling.
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
    U = lu.matrixLU().triangularView<Eigen::Upper>();
    P = lu.permutationP().inverse();
    Q = lu.permutationQ().inverse();
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
    Vector x = P * b;
    x = L.topLeftCorner(rank, rank).triangularView<Eigen::Lower>().solve(x);
    x = U.topRows(rank).triangularView<Eigen::Upper>().solve(x);
    x = Q * x;
    return x;
}

auto LU::trsolve(const Vector& b) -> Vector
{
    const auto& indices = Q.indices();
    Vector x(rank);
    for(Index i = 0; i < rank; ++i)
        x[i] = b[indices[i]];
    x = tr(U).topLeftCorner(rank, rank).triangularView<Eigen::Lower>().solve(x);
    x = tr(L).topRows(rank).triangularView<Eigen::Upper>().solve(x);
    x = P.inverse() * x;
    return x;
}

} // namespace Reaktoro
