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
#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/LU>
using namespace Eigen;

namespace Reaktor {
namespace {

auto toeigen(const Matrix& mat) -> MatrixXd
{
    MatrixXd res(mat.n_rows, mat.n_cols);
    for(unsigned i = 0; i < res.rows(); ++i)
        for(unsigned j = 0; j < res.cols(); ++j)
            res(i, j) = mat(i, j);
    return res;
}

auto toeigen(const Vector& vec) -> VectorXd
{
    VectorXd res(vec.n_rows);
    for(unsigned i = 0; i < res.rows(); ++i)
        res[i] = vec[i];
    return res;
}

auto toarma(const MatrixXd& mat) -> Matrix
{
    Matrix res(mat.rows(), mat.cols());
    for(unsigned i = 0; i < res.n_rows; ++i)
        for(unsigned j = 0; j < res.n_cols; ++j)
            res(i, j) = mat(i, j);
    return res;
}

auto toeigen(const VectorXd& vec) -> Vector
{
    Vector res(vec.rows());
    for(unsigned i = 0; i < res.n_rows; ++i)
        res[i] = vec[i];
    return res;
}

} // namespace

auto solveFull(const SaddlePointProblem& problem, SaddlePointResult& result) -> void
{
    const MatrixXd A = toeigen(problem.A);
    const MatrixXd H = toeigen(problem.H);
    const MatrixXd f = toeigen(problem.f);
    const MatrixXd g = toeigen(problem.g);

    const unsigned m = A.rows();
    const unsigned n = A.cols();

    MatrixXd lhs(n + m, n + m);
    lhs.block(0, 0, n, n).noalias() = H;
    lhs.block(0, n, n, m).noalias() = -A.transpose();
    lhs.block(n, 0, m, n).noalias() = A;
    lhs.block(n, n, m, m).noalias() = MatrixXd::Zero(m, m);

    VectorXd rhs(n + m);
    rhs.segment(0, n) = f;
    rhs.segment(n, m) = g;

    VectorXd sol = lhs.fullPivLu().solve(rhs);
    VectorXd x = sol.segment(0, n);
    VectorXd y = sol.segment(n, m);

    assert((lhs*sol-rhs).norm()/rhs.norm() < 1e-13);

    result.solution.x = toarma(x);
    result.solution.y = toarma(y);
}

auto solveNullspace(const SaddlePointProblem& problem, SaddlePointResult& result) -> void
{
    const MatrixXd A = toeigen(problem.A);
    const MatrixXd H = toeigen(problem.H);
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
