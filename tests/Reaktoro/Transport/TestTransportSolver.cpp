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

// C++ includes
#include <iostream>

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

TEST_CASE("Testing Mesh")
{
    Mesh mesh;
    mesh.setDiscretization(10, 1.0, 2.0);

    CHECK(mesh.numCells() == 10);
    CHECK(mesh.xl() == 1.0);
    CHECK(mesh.xr() == 2.0);
    CHECK(mesh.dx() == approx(0.1));
    CHECK(mesh.xcells().isApprox(linspace(10, 1.05, 1.95)));
}

TEST_CASE("Testing transport solver for a pure advection problem")
{
    const auto num_cells = 10;
    const auto num_steps = 100;
    const auto velocity = 1.0;

    Mesh mesh(num_cells);

    const auto dx = mesh.dx();
    const auto xcells = mesh.xcells();
    const double cfl_values[] = { 0.1, 0.5, 1.0 };

    for(double cfl : cfl_values)
    {
        TransportSolver transport;
        transport.setMesh(mesh);
        transport.setBoundaryCondition(1.0);
        transport.setVelocity(velocity);

        Vector u = zeros(num_cells);
        Vector expected = zeros(num_cells);

        double t = 0.0;
        const double dt = cfl*dx/velocity;

        transport.setTimeStep(dt);

        for(Index i = 0; i < num_steps; ++i)
        {
            transport.step(u);
            t += dt;
            expected.head(std::min(int(t/dx), num_cells)).fill(1.0);

//            std::cout << "u(actual)   = " << tr(u) << std::endl;
//            std::cout << "u(expected) = " << tr(expected) << std::endl;
//            std::cout << "norm(error) = " << norm(u - expected) << std::endl;

            CHECK(norminf(u - expected) <= 1.0);
        }
    }
}

TEST_CASE("Testing transport solver")
{
    const auto num_cells = 100;
    const auto num_steps = 20;
    const auto diffusion = 1.0e-5;
    const auto velocity = 1.0;

    Mesh mesh(num_cells);

    const auto dx = mesh.dx();
    const auto xcells = mesh.xcells();
    const auto sinx = Vector(xcells.array().sin());
    const auto cosx = Vector(xcells.array().cos());

    TransportSolver transport;
    transport.setMesh(mesh);
    transport.setBoundaryCondition(0.0);
    transport.setDiffusionCoeff(diffusion);
    transport.setVelocity(velocity);

    auto update_source = [=](VectorRef q, double t)
    {
        q.noalias() = (1 + diffusion*t)*sinx + velocity*t*cosx;
    };

    Vector q = zeros(num_cells);
    Vector u = zeros(num_cells);
    Vector expected = zeros(num_cells);
    double t = 0.0;

    SUBCASE("When CFL = v*dt/dx = 1.0")
    {
        const double cfl = 0.1;
        const double dt = cfl*dx/velocity;

        transport.setTimeStep(dt);

        for(Index i = 0; i < num_steps; ++i)
        {
            update_source(q, t);
            transport.step(u, q);
            t += dt;
            expected.noalias() = t*sinx;

//            std::cout << "u(actual)   = " << tr(u) << std::endl;
//            std::cout << "u(expected) = " << tr(expected) << std::endl;
//            std::cout << "norm(error) = " << norm(u - expected) << std::endl;

            CHECK(norminf(u - expected) < 1e-2);
        }
    }
}
