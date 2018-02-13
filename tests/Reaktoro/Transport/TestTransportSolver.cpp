// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2017 Allan Leal
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

#include <doctest/doctest.hpp>

// Reaktoro includes
#include <Reaktoro/Math/Eigen/LU>
#include <Reaktoro/Transport/TransportSolver.hpp>
using namespace Reaktoro;

TEST_CASE("Testing tridiagonal solver")
{
    const Index n = 10;

    TridiagonalMatrix trimat(n);
    trimat.a() = abs(random(n - 1));
    trimat.b() = 2 * abs(random(n));
    trimat.c() = abs(random(n - 1));

    Matrix mat = zeros(n, n);
    mat.diagonal(-1) = trimat.a();
    mat.diagonal( 0) = trimat.b();
    mat.diagonal(+1) = trimat.c();

    CHECK(mat.isApprox(Matrix(trimat)));

    Vector expected = linspace(n, 0, n - 1);

    Vector d = mat * expected;

    Vector x(n);

    trimat.factorize();
    trimat.solve(x, d);

    CHECK(x.isApprox(expected));

    x = d;
    trimat.solve(x);

    CHECK(x.isApprox(expected));
}
