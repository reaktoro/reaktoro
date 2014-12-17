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

// Eigen includes
#include <Reaktor/eigen/Cholesky>
#include <Reaktor/eigen/Core>
#include <Reaktor/eigen/LU>
using namespace Eigen;

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Common/TimeUtils.hpp>

namespace Reaktor {

auto solve(const SaddlePointProblem& problem, SaddlePointResult& result, const SaddlePointOptions& options) -> void
{
    switch(options.algorithm)
    {
    case FullspaceDense    : solveFullspaceDense(problem, result, options); break;
    case FullspaceSparse   : solveFullspaceSparse(problem, result, options); break;
    case Rangespace        : solveRangespace(problem, result, options); break;
    case Nullspace         : solveNullspace(problem, result, options); break;
    case NullspacePartial  : solveNullspacePartial(problem, result, options); break;
    default: error("Cannot solve the saddle point problem.", "Unknown algorithm.");
    }
}

auto solveFullspaceDense(const SaddlePointProblem& problem, SaddlePointResult& result, const SaddlePointOptions& options) -> void
{
    Time begin = time();

    const auto& A = problem.A;
    const auto& H = problem.H;
    const auto& f = problem.f;
    const auto& g = problem.g;

    auto& lhs = result.internal.lhs;
    auto& rhs = result.internal.rhs;
    auto& sol = result.internal.sol;

    const unsigned m = A.rows();
    const unsigned n = A.cols();

    lhs.resize(n + m, n + m);
    rhs.resize(n + m);
    sol.resize(n + m);

    lhs.block(0, 0, n, n).noalias() =  H;
    lhs.block(0, n, n, m).noalias() = -A.transpose();
    lhs.block(n, 0, m, n).noalias() =  A;
    lhs.block(n, n, m, m).noalias() =  zeros(m, m);

    rhs.segment(0, n) = f;
    rhs.segment(n, m) = g;

    // Solve with full partial pivoting for higher stability
    sol = lhs.fullPivLu().solve(rhs);

    result.solution.x = rows(sol, 0, n);
    result.solution.y = rows(sol, n, n + m);

    Time end = time();
    result.statistics.converged = true;
    result.statistics.time = elapsed(end, begin);
}

auto solveFullspaceSparse(const SaddlePointProblem& problem, SaddlePointResult& result, const SaddlePointOptions& options) -> void
{
    Time begin = time();

    error("Cannot solve the saddle point problem.", "The FullspaceSparse algorithm has not been implemented yet.");

    Time end = time();
    result.statistics.converged = true;
    result.statistics.time = elapsed(end, begin);
}

auto solveRangespace(const SaddlePointProblem& problem, SaddlePointResult& result, const SaddlePointOptions& options) -> void
{
    Time begin = time();

    const auto& H = problem.H;
    const auto& A = problem.A;
    const auto& f = problem.f;
    const auto& g = problem.g;
    auto& x = result.solution.x;
    auto& y = result.solution.y;

    if(options.properties & DiagonalH)
    {
        const Vector invH = H.diagonal().cwiseInverse();
        const Matrix AinvH = A * invH.asDiagonal();
        const Matrix AinvHAT = AinvH * A.transpose();
        y = AinvHAT.llt().solve(g - AinvH*f);
        x = invH.asDiagonal() * f + AinvH.transpose()*y;
    }
    else
    {
        error("Cannot solve the saddle point problem.", "The Rangespace algorithm has been implemented only for diagonal matrices H.");
    }

    Time end = time();
    result.statistics.converged = true;
    result.statistics.time = elapsed(end, begin);
}

auto solveNullspace(const SaddlePointProblem& problem, SaddlePointResult& result, const SaddlePointOptions& options) -> void
{
    Time begin = time();

    const auto& A = problem.A;
    const auto& H = problem.H;
    const auto& f = problem.f;
    const auto& g = problem.g;

    const unsigned m = A.rows();
    const unsigned n = A.cols();

    const Matrix Heff = 0.5 * (H + H.transpose());

    const FullPivLU<Matrix> lu_A = A.fullPivLu();
    const Matrix L = lu_A.matrixLU().leftCols(m).triangularView<UnitLower>();
    const Matrix U = lu_A.matrixLU().triangularView<Upper>();
    const auto P1 = lu_A.permutationP();
    const auto P2 = lu_A.permutationQ();
    const auto U1 = U.leftCols(m);
    const auto U2 = U.rightCols(n - m);

    Matrix Z(n, n - m);
    Z.topRows(m) = -U1.triangularView<Upper>().solve(U2);
    Z.bottomRows(n - m) = Matrix::Identity(n - m, n - m);

    Z = P2*Z;

    Matrix Y = Matrix::Zero(n, m);
    Y.topRows(m) = L.triangularView<Lower>().solve(Matrix::Identity(m, m));
    Y.topRows(m) = U1.triangularView<Upper>().solve(Y.topRows(m));

    Y = P2*Y*P1;

    Matrix ZTHZ = Z.transpose() * Heff * Z;

    LLT<Matrix> llt_ZTHZ(ZTHZ);

    Vector tmp = Z.transpose() * (f - Heff*Y*g);

    Vector xZ = llt_ZTHZ.solve(tmp);

    result.solution.x = Z*xZ + Y*g;
    result.solution.y = Y.transpose() * (Heff*result.solution.x - f);

    Time end = time();
    result.statistics.converged = true;
    result.statistics.time = elapsed(end, begin);
}

auto solveNullspacePartial(const SaddlePointProblem& problem, SaddlePointResult& result, const SaddlePointOptions& options) -> void
{
    Time begin = time();

    error("Cannot solve the saddle point problem.", "The NullspacePartial algorithm has not been implemented yet.");

    Time end = time();
    result.statistics.converged = true;
    result.statistics.time = elapsed(end, begin);
}

} // namespace Reaktor
