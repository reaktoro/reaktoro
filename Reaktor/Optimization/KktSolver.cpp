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

#include "KktSolver.hpp"

// Eigen includes
#include <Reaktor/eigen/Cholesky>
#include <Reaktor/eigen/Core>
#include <Reaktor/eigen/LU>
using namespace Eigen;

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Common/TimeUtils.hpp>

namespace Reaktor {

struct KktProblemDense
{
    /// The internal data for the KKT problem
    Matrix kkt_lhs;
    Vector kkt_rhs;
    Vector kkt_sol;
    PartialPivLU<Matrix> kkt_lu;

    /// Decompose any necessary matrix before the KKT calculation.
    /// Note that this method should be called before `solve`,
    /// once the matrices `H` and `A` have been initialised.
    /// @see solve
    auto decompose(const KktMatrix& lhs) -> void;

    /// Solve the KKT problem using a dense LU decomposition.
    /// Note that this method requires `decompose` to be called a priori.
    /// @param rhs The right-hand side vector of the KKT equation
    /// @param result[out] The result of the KKT equation
    /// @see decompose
    auto solve(const KktMatrix& lhs, const KktVector& rhs, KktResult& result) -> void;
};

struct KktProblemInverseH
{
    /// The internal data for the KKT problem
    Matrix AinvH;
    Matrix AinvHAt;
    LLT<Matrix> llt_AinvHAt;

    /// Decompose any necessary matrix before the KKT calculation.
    /// Note that this method should be called before `solve`,
    /// once the matrices `H` and `A` have been initialised.
    /// @see solve
    auto decompose(const KktMatrix& lhs) -> void;

    /// Solve the KKT problem using an efficient rangespace decomposition approach.
    /// Note that this method requires `decompose` to be called a priori.
    /// @param rhs The right-hand side vector of the KKT equation
    /// @param result[out] The result of the KKT equation
    /// @see decompose
    auto solve(const KktMatrix& lhs, const KktVector& rhs, KktResult& result) -> void;
};

struct KktProblemDiagonalH
{
    /// The internal data for the KKT problem
    Matrix AinvH;
    Matrix AinvHAt;
    LLT<Matrix> llt_AinvHAt;

    /// Decompose any necessary matrix before the KKT calculation.
    /// Note that this method should be called before `solve`,
    /// once the matrices `H` and `A` have been initialised.
    /// @see solve
    auto decompose(const KktMatrix& lhs) -> void;

    /// Solve the KKT problem using an efficient rangespace decomposition approach.
    /// Note that this method requires `decompose` to be called a priori.
    /// @param rhs The right-hand side vector of the KKT equation
    /// @param result[out] The result of the KKT equation
    /// @see decompose
    auto solve(const KktMatrix& lhs, const KktVector& rhs, KktResult& result) -> void;
};

struct KktProblemConstantA
{
    /// The nullspace matrix of `A` with the property `AZ = 0`
    Matrix Z;

    /// The rangespace matrix of `A` with the property `AY = I`
    Matrix Y;

    /// Auxiliary data for the nullspace algorithm
    Matrix ZtHZ;
    LLT<Matrix> llt_ZtHZ;
    Vector xZ;

    /// Auxiliary data for finding the nullspace and rangespace matrices `Z` and `Y`
    FullPivLU<Matrix> lu_A;
    Matrix L;
    Matrix U;

    /// Initialise the constant bottom-left matrix `A` of the KKT equation.
    /// This method should be called once to initialise the `A` matrix
    /// of the KKT equation and whenever it is changed subsequently.
    auto initialise(const Matrix& A) -> void;

    /// Decompose any necessary matrix before the KKT calculation.
    /// Note that this method should be called before `solve`,
    /// once the matrices `H` and `A` have been initialised.
    /// @see solve
    auto decompose(const KktMatrix& lhs) -> void;

    /// Solve the KKT problem using an efficient nullspace decomposition approach.
    /// Note that this method requires `decompose` to be called a priori.
    /// @param rhs The right-hand side vector of the KKT equation
    /// @param result[out] The result of the KKT equation
    /// @see decompose
    auto solve(const KktMatrix& lhs, const KktVector& rhs, KktResult& result) -> void;
};

auto KktProblemDense::decompose(const KktMatrix& lhs) -> void
{
    // The references to the H and A matrices
    const auto& H = lhs.H;
    const auto& A = lhs.A;

    // The dimensions of the KKT problem
    const unsigned n = A.cols();
    const unsigned m = A.rows();

    // Ensure the components of the KKT equation have adequate dimensions
    kkt_lhs.resize(n + m, n + m);
    kkt_rhs.resize(n + m);
    kkt_sol.resize(n + m);

    // Assemble the left-hand side of the KKT equation
    kkt_lhs.block(0, 0, n, n).noalias() =  H;
    kkt_lhs.block(0, n, n, m).noalias() = -A.transpose();
    kkt_lhs.block(n, 0, m, n).noalias() =  A;
    kkt_lhs.block(n, n, m, m).noalias() =  zeros(m, m);

    // Perform the LU decomposition
    kkt_lu.compute(kkt_lhs);
}

auto KktProblemDense::solve(const KktMatrix& lhs, const KktVector& rhs, KktResult& result) -> void
{
    // The dimensions of the KKT problem
    const unsigned n = lhs.A.cols();
    const unsigned m = lhs.A.rows();

    // Assemble the right-hand side of the KKT equation
    kkt_rhs.segment(0, n) = rhs.f;
    kkt_rhs.segment(n, m) = rhs.g;

    // Check if the LU decomposition has already been performed
    if(kkt_lu.rows() != n + m or kkt_lu.cols() != n + m)
        error("Cannot solve the KKT equation using the fullspace algorithm.",
            "The LU decomposition of the KKT matrix has not been initialised"
            "or updated for a new problem with different dimension.");

    // Solve the linear system with the LU decomposition already calculated
    kkt_sol = kkt_lu.solve(kkt_rhs);

    // Extract the solution `x` and `y` from the linear system solution `sol`
    result.solution.x = rows(kkt_sol, 0, n);
    result.solution.y = rows(kkt_sol, n, m);

    // Set the statistics of the calculation
    result.statistics.converged = true;
}

auto KktProblemInverseH::decompose(const KktMatrix& lhs) -> void
{
    AinvH = lhs.A * lhs.invH;
    AinvHAt = AinvH * lhs.A.transpose();

    llt_AinvHAt.compute(AinvHAt);

    if(llt_AinvHAt.info() != Eigen::Success)
        error("Cannot solve the KKT problem with the rangespace algorithm.",
            "The provided matrix `H` of the KKT equation might not be symmetric positive-definite.");
}

auto KktProblemInverseH::solve(const KktMatrix& lhs, const KktVector& rhs, KktResult& result) -> void
{
    // The references to the right-hand side components of the KKT problem
    const auto& f = rhs.f;
    const auto& g = rhs.g;

    // The references to the solution of the KKT problem
    auto& x = result.solution.x;
    auto& y = result.solution.y;

    y = llt_AinvHAt.solve(g - AinvH*f);
    x = lhs.invH * f + AinvH.transpose()*y;

    result.statistics.converged = true;
}

auto KktProblemDiagonalH::decompose(const KktMatrix& lhs) -> void
{
    const auto invH = inv(lhs.diagH);
    AinvH = lhs.A * diag(invH);
    AinvHAt = AinvH * lhs.A.transpose();

    llt_AinvHAt.compute(AinvHAt);

    if(llt_AinvHAt.info() != Eigen::Success)
        error("Cannot solve the KKT problem with the rangespace algorithm.",
            "The provided matrix `H` of the KKT equation might not be symmetric positive-definite.");
}

auto KktProblemDiagonalH::solve(const KktMatrix& lhs, const KktVector& rhs, KktResult& result) -> void
{
    // The references to the right-hand side components of the KKT problem
    const auto& f = rhs.f;
    const auto& g = rhs.g;

    // The references to the solution of the KKT problem
    auto& x = result.solution.x;
    auto& y = result.solution.y;

    y = llt_AinvHAt.solve(g - AinvH*f);
    x = inv(lhs.diagH) % f + AinvH.transpose()*y;

    result.statistics.converged = true;
}

auto KktProblemConstantA::initialise(const Matrix& A) -> void
{
    // The dimensions of the matrix `A`
    const unsigned n = A.cols();
    const unsigned m = A.rows();

    // Perform a LU decomposition of the matrix `A`
    lu_A.compute(A);

    // Get the lower and upper matrices
    L = lu_A.matrixLU().leftCols(m).triangularView<UnitLower>();
    U = lu_A.matrixLU().triangularView<Upper>();

    // Get the permutation matrices
    const auto P1 = lu_A.permutationP();
    const auto P2 = lu_A.permutationQ();

    // Set the U1 and U2 submatrices of U = [U1 U2]
    const auto U1 = U.leftCols(m);
    const auto U2 = U.rightCols(n - m);

    // Update the nullspace matrix `Z` of `A`
    Z = zeros(n, n - m);
    Z.topRows(m) = -U1.triangularView<Upper>().solve(U2);
    Z.bottomRows(n - m) = identity(n - m, n - m);
    Z = P2*Z;

    // Update the rangespace matrix `Y` of `A`
    Y = zeros(n, m);
    Y.topRows(m) = L.triangularView<Lower>().solve(identity(m, m));
    Y.topRows(m) = U1.triangularView<Upper>().solve(Y.topRows(m));
    Y = P2*Y*P1;
}

auto KktProblemConstantA::decompose(const KktMatrix& lhs) -> void
{
    // Check if the nullspace and rangespace matrices have been initialised
    if(not Z.size())
        error("Cannot solve the KKT equation using the nullspace algorithm.",
            "The matrix `A` of the KKT equation was not initialised.");

    // Compute the reduced nullspace matrix
    ZtHZ = Z.transpose() * lhs.H * Z;

    // Compute Cholesky decomposition of the reduced nullspace matrix
    llt_ZtHZ.compute(ZtHZ);

    // Check if the Cholesky decomposition was successful
    if(llt_ZtHZ.info() != Eigen::Success)
        error("Cannot solve the KKT equation using the nullspace algorithm.",
            "The provided H matrix might not be symmetric positive-definite or is ill-conditioned.");
}

auto KktProblemConstantA::solve(const KktMatrix& lhs, const KktVector& rhs, KktResult& result) -> void
{
    // The references to the right-hand side components of the KKT problem
    const auto& f = rhs.f;
    const auto& g = rhs.g;

    // The references to the solution of the KKT problem
    auto& x = result.solution.x;
    auto& y = result.solution.y;

    // The dimensions of `x` and `y`
    const unsigned n = x.rows();
    const unsigned m = y.rows();

    // Check if the Cholesky decomposition has already been performed
    if(llt_ZtHZ.rows() != n - m or llt_ZtHZ.cols() != n - m)
        error("Cannot solve the KKT equation using the nullspace algorithm.",
            "The Cholesky decomposition of the reduced Hessian matrix has not "
            "been initialised or updated for a new problem with different dimension.");

    // Compute the `xZ` component of `x`
    xZ = Z.transpose() * (f - lhs.H*Y*g);
    xZ = llt_ZtHZ.solve(xZ);

    // Compute both `x` and `y` variables
    x = Z*xZ + Y*g;
    y = Y.transpose() * (lhs.H*x - f);

    result.statistics.converged = true;
}

struct KktProblem::Impl
{
    KktProblemDense kkt_dense;
    KktProblemConstantA kkt_const_A;
    KktProblemDiagonalH kkt_diag_H;
    KktProblemInverseH kkt_inv_H;
};

KktProblem::KktProblem()
: pimpl(new Impl())
{}

KktProblem::KktProblem(const KktProblem& other)
: pimpl(new Impl(*other.pimpl))
{}

KktProblem::~KktProblem()
{}

auto KktProblem::operator=(KktProblem other) -> KktProblem&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto KktProblem::setConstantA(const Matrix& A) -> void
{
    pimpl->kkt_const_A.initialise(A);
}

auto KktProblem::decompose(const KktMatrix& lhs) -> void
{
    pimpl->kkt_dense.decompose(lhs);
}

auto KktProblem::decomposeWithInverseH(const KktMatrix& lhs) -> void
{
    pimpl->kkt_inv_H.decompose(lhs);
}

auto KktProblem::decomposeWithDiagonalH(const KktMatrix& lhs) -> void
{
    pimpl->kkt_diag_H.decompose(lhs);
}

auto KktProblem::decomposeWithConstantA(const KktMatrix& lhs) -> void
{
    pimpl->kkt_const_A.decompose(lhs);
}

auto KktProblem::solve(const KktMatrix& lhs, const KktVector& rhs, KktResult& result) -> void
{
    pimpl->kkt_dense.solve(lhs, rhs, result);
}

auto KktProblem::solveWithInverseH(const KktMatrix& lhs, const KktVector& rhs, KktResult& result) -> void
{
    pimpl->kkt_inv_H.solve(lhs, rhs, result);
}

auto KktProblem::solveWithDiagonalH(const KktMatrix& lhs, const KktVector& rhs, KktResult& result) -> void
{
    pimpl->kkt_diag_H.solve(lhs, rhs, result);
}

auto KktProblem::solveWithConstantA(const KktMatrix& lhs, const KktVector& rhs, KktResult& result) -> void
{
    pimpl->kkt_const_A.solve(lhs, rhs, result);
}

} // namespace Reaktor
