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

#include "SaddlePointUtils.hpp"

namespace Reaktor {

auto solveNullspace(const SaddlePointProblem& problem, SaddlePointResult& result) -> void
{
    const unsigned m = problem.A.n_rows;
    const unsigned n = problem.A.n_cols;

    const Matrix& H = problem.H;
    const Matrix& A = problem.A;
    const Vector& f = problem.f;
    const Vector& g = problem.g;

    Matrix& Z    = result.internal.Z;
    Matrix& Y    = result.internal.Y;
    Matrix& ZtHZ = result.internal.ZtHZ;
    Vector& xZ   = result.internal.xZ;
    Matrix& L    = result.internal.L;
    Matrix& U    = result.internal.U;
    Matrix& P    = result.internal.P;
    Matrix& R    = result.internal.R;
    Vector& x    = result.solution.x;
    Vector& y    = result.solution.y;

    // Compute the LU decomposition of matrix A as A = P.t() * L * U
    arma::lu(L, U, P, A);

    // Extract the U1 and U2 part of U, where U = [U1 U2]
    const auto U1 = U.cols(0, m - 1);
    const auto U2 = U.cols(m, n - 1);

    // Compute the Z matrix, where AZ = 0
    Z.resize(n, n - m);
    Z.rows(0, m - 1) = arma::solve(arma::trimatu(-U1), U2);
    Z.rows(m, n - 1) = arma::eye(n - m, n - m);

    // Compute the Y matrix, where AY = I
    Y = zeros(n, m);
    Y.rows(0, m - 1) = arma::solve(arma::trimatl(L), P);
    Y.rows(0, m - 1) = arma::solve(arma::trimatu(U1), Y.rows(0, m - 1));

    ZtHZ = Z.t() * H * Z;

    R = chol(ZtHZ);
    xZ = Z.t() * (f - H*Y*g);
    arma::solve(xZ, arma::trimatl(R.t()), xZ);
    arma::solve(xZ, arma::trimatu(R), xZ);

    x = Z*xZ + Y*g;
    y = Y.t() * (H*x - f);
}

} // namespace Reaktor
