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
namespace {

auto toeigen(const Matrix& mat) -> Eigen::Map<const Eigen::MatrixXd>
{
    return Eigen::Map<const Eigen::MatrixXd>(mat.memptr(), mat.n_rows, mat.n_cols);
}

auto toeigen(const Vector& vec) -> Eigen::Map<const Eigen::VectorXd>
{
    return Eigen::Map<const Eigen::VectorXd>(vec.memptr(), vec.n_rows);
}

auto toeigen(Matrix& mat) -> Eigen::Map<Eigen::MatrixXd>
{
    return Eigen::Map<Eigen::MatrixXd>(mat.memptr(), mat.n_rows, mat.n_cols);
}

auto toeigen(Vector& vec) -> Eigen::Map<Eigen::VectorXd>
{
    return Eigen::Map<Eigen::VectorXd>(vec.memptr(), vec.n_rows);
}

auto toarma(const MatrixXd& mat) -> Matrix
{
    Matrix res(mat.rows(), mat.cols());
    for(unsigned i = 0; i < res.n_rows; ++i)
        for(unsigned j = 0; j < res.n_cols; ++j)
            res(i, j) = mat(i, j);
    return res;
}

auto toarma(const VectorXd& vec) -> Vector
{
    Vector res(vec.rows());
    for(unsigned i = 0; i < res.n_rows; ++i)
        res[i] = vec[i];
    return res;
}

} // namespace

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

    const auto A = toeigen(problem.A);
    const auto H = toeigen(problem.H);
    const auto f = toeigen(problem.f);
    const auto g = toeigen(problem.g);

    const unsigned m = A.rows();
    const unsigned n = A.cols();

    result.internal.lhs.resize(n + m, n + m);
    result.internal.rhs.resize(n + m);
    result.internal.sol.resize(n + m);

    auto lhs = toeigen(result.internal.lhs);
    auto rhs = toeigen(result.internal.rhs);
    auto sol = toeigen(result.internal.sol);

    lhs.block(0, 0, n, n).noalias() = H;
    lhs.block(0, n, n, m).noalias() = -A.transpose();
    lhs.block(n, 0, m, n).noalias() = A;
    lhs.block(n, n, m, m).noalias() = MatrixXd::Zero(m, m);

    rhs.segment(0, n) = f;
    rhs.segment(n, m) = g;

    // Solve with full partial pivoting for higher stability
    sol = lhs.fullPivLu().solve(rhs);

    result.solution.x = result.internal.sol.subvec(0, n-1);
    result.solution.y = result.internal.sol.subvec(n, n+m-1);

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
        const Vector invH = 1/H.diag();
        const Matrix AinvH = A * arma::diagmat(invH);
        const Matrix AinvHAT = AinvH * A.t();
        Matrix R = arma::chol(AinvHAT);
        y = g - AinvH*f;
        arma::solve(y, arma::trimatl(R.t()), y);
        arma::solve(y, arma::trimatu(R), y);
        x = invH % f + AinvH.t()*y;
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

    const auto A = toeigen(problem.A);
    const auto H = toeigen(problem.H);
    const MatrixXd f = toeigen(problem.f);
    const MatrixXd g = toeigen(problem.g);

    const unsigned m = A.rows();
    const unsigned n = A.cols();

    const MatrixXd Heff = 0.5 * (H + H.transpose());

    const FullPivLU<MatrixXd> lu_A = A.fullPivLu();
    const MatrixXd L = lu_A.matrixLU().leftCols(m).triangularView<UnitLower>();
    const MatrixXd U = lu_A.matrixLU().triangularView<Upper>();
    const auto P1 = lu_A.permutationP();
    const auto P2 = lu_A.permutationQ();
    const auto U1 = U.leftCols(m);
    const auto U2 = U.rightCols(n - m);

    MatrixXd Z(n, n - m);
    Z.topRows(m) = -U1.triangularView<Upper>().solve(U2);
    Z.bottomRows(n - m) = MatrixXd::Identity(n - m, n - m);

    Z = P2*Z;

    MatrixXd Y = MatrixXd::Zero(n, m);
    Y.topRows(m) = L.triangularView<Lower>().solve(MatrixXd::Identity(m, m));
    Y.topRows(m) = U1.triangularView<Upper>().solve(Y.topRows(m));

    Y = P2*Y*P1;

    MatrixXd ZTHZ = Z.transpose() * Heff * Z;

    LLT<MatrixXd> llt_ZTHZ(ZTHZ);

    VectorXd tmp = Z.transpose() * (f - Heff*Y*g);

    VectorXd xZ = llt_ZTHZ.solve(tmp);

    VectorXd x = Z*xZ + Y*g;
    VectorXd y = Y.transpose() * (Heff*x - f);

    result.solution.x = toarma(x);
    result.solution.y = toarma(y);

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

//// Reaktor is a C++ library for computational reaction modelling.
////
//// Copyright (C) 2014 Allan Leal
////
//// This program is free software: you can redistribute it and/or modify
//// it under the terms of the GNU General Public License as published by
//// the Free Software Foundation, either version 3 of the License, or
//// (at your option) any later version.
////
//// This program is distributed in the hope that it will be useful,
//// but WITHOUT ANY WARRANTY; without even the implied warranty of
//// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//// GNU General Public License for more details.
////
//// You should have received a copy of the GNU General Public License
//// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
//#include "SaddlePointUtils.hpp"
//
//namespace Reaktor {
//
//auto solveNullspace(const SaddlePointProblem& problem, SaddlePointResult& result) -> void
//{
//    const unsigned m = problem.A.n_rows;
//    const unsigned n = problem.A.n_cols;
//
//    const Matrix& H = problem.H;
//    const Matrix& A = problem.A;
//    const Vector& f = problem.f;
//    const Vector& g = problem.g;
//
//    Matrix& Z    = result.internal.Z;
//    Matrix& Y    = result.internal.Y;
//    Matrix& ZtHZ = result.internal.ZtHZ;
//    Vector& xZ   = result.internal.xZ;
//    Matrix& L    = result.internal.L;
//    Matrix& U    = result.internal.U;
//    Matrix& P    = result.internal.P;
//    Matrix& R    = result.internal.R;
//    Vector& x    = result.solution.x;
//    Vector& y    = result.solution.y;
//
//    // Compute the LU decomposition of matrix A as A = P.t() * L * U
//    arma::lu(L, U, P, A);
//
//    // Extract the U1 and U2 part of U, where U = [U1 U2]
//    const auto U1 = U.cols(0, m - 1);
//    const auto U2 = U.cols(m, n - 1);
//
//    // Compute the Z matrix, where AZ = 0
//    Z.resize(n, n - m);
//    Z.rows(0, m - 1) = arma::solve(arma::trimatu(-U1), U2);
//    Z.rows(m, n - 1) = arma::eye(n - m, n - m);
//
//    // Compute the Y matrix, where AY = I
//    Y = zeros(n, m);
//    Y.rows(0, m - 1) = arma::solve(arma::trimatl(L), P);
//    Y.rows(0, m - 1) = arma::solve(arma::trimatu(U1), Y.rows(0, m - 1));
//
//    ZtHZ = Z.t() * H * Z;
//
//    R = chol(ZtHZ);
//    xZ = Z.t() * (f - H*Y*g);
//    arma::solve(xZ, arma::trimatl(R.t()), xZ);
//    arma::solve(xZ, arma::trimatu(R), xZ);
//
//    x = Z*xZ + Y*g;
//    y = Y.t() * (H*x - f);
//}
//
//} // namespace Reaktor
