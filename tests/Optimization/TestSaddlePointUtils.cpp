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

#include "TestSaddlePointUtils.hpp"

// Reaktoro includes
#include <Reaktoro/Reaktoro.hpp>

namespace Reaktoro {
namespace {

auto test_solveFullspaceDense() -> void
{
    const unsigned m = 3;
    const unsigned n = 7;

    const Matrix M = arma::randu(n, n);

    SaddlePointProblem problem;
    problem.H = M*M.t();
    problem.A = arma::randu(m, n);
    problem.f = arma::randu(n);
    problem.g = arma::randu(m);

    SaddlePointResult state;
    SaddlePointOptions options;

    solveFullspaceDense(problem, state, options);

    const auto& H = problem.H;
    const auto& A = problem.A;
    const auto& f = problem.f;
    const auto& g = problem.g;
    const auto& x = state.x;
    const auto& y = state.y;
    const double residualx = arma::norm(H*x - A.t()*y - f)/arma::norm(f);
    const double residualy = arma::norm(A*x - g)/arma::norm(g);

    ASSERT(residualx < 1e-14);
    ASSERT(residualy < 1e-14);
}

auto test_solveNullspace() -> void
{
    const unsigned m = 3;
    const unsigned n = 7;

    const Matrix M = arma::randu(n, n);

    SaddlePointProblem problem;
    problem.H = M*M.t();
    problem.A = arma::randu(m, n);
    problem.f = arma::randu(n);
    problem.g = arma::randu(m);

    SaddlePointResult state;
    SaddlePointOptions options;

    solveNullspace(problem, state, options);

    const auto& H = problem.H;
    const auto& A = problem.A;
    const auto& f = problem.f;
    const auto& g = problem.g;
    const auto& x = state.x;
    const auto& y = state.y;
    const double residualx = arma::norm(H*x - A.t()*y - f)/arma::norm(f);
    const double residualy = arma::norm(A*x - g)/arma::norm(g);

    ASSERT(residualx < 1e-14);
    ASSERT(residualy < 1e-14);
}

auto test_solveRangespaceDiagonal() -> void
{
    const unsigned m = 3;
    const unsigned n = 7;

    const Matrix M = arma::randu(n, n) + 1.0;

    SaddlePointProblem problem;
    problem.H = arma::diagmat(M);
    problem.A = arma::randu(m, n);
    problem.f = arma::randu(n);
    problem.g = arma::randu(m);

    SaddlePointResult state;
    SaddlePointOptions options;
    options.properties = DiagonalH;

    solveRangespace(problem, state, options);

    const auto& H = problem.H;
    const auto& A = problem.A;
    const auto& f = problem.f;
    const auto& g = problem.g;
    const auto& x = state.x;
    const auto& y = state.y;
    const double residualx = arma::norm(H*x - A.t()*y - f)/arma::norm(f);
    const double residualy = arma::norm(A*x - g)/arma::norm(g);

    ASSERT(residualx < 1e-14);
    ASSERT(residualy < 1e-14);
}

} // namespace

auto testSuiteSaddlePointUtils() -> cute::suite
{
    cute::suite s;

    s += CUTE(test_solveFullspaceDense);
    s += CUTE(test_solveNullspace);
    s += CUTE(test_solveRangespaceDiagonal);

    return s;
}

} // namespace Reaktoro
