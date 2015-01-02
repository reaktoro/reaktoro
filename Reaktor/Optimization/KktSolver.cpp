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

struct KktProblemFull
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
    auto decompose(const Matrix& H, const Matrix& A) -> void;

    /// Solve the KKT problem using a dense LU decomposition.
    /// Note that this method requires `decompose` to be called a priori.
    /// @param rhs The right-hand side vector of the KKT equation
    /// @param result[out] The result of the KKT equation
    /// @see decompose
    auto solve(const Vector& a, const Vector& b, Vector& x, Vector& y) -> void;
};

struct KktProblemInverseH
{
    /// The internal data for the KKT problem
    const Matrix* invH_ptr;
    Matrix AinvH;
    Matrix AinvHAt;
    LLT<Matrix> llt_AinvHAt;

    /// Decompose any necessary matrix before the KKT calculation.
    /// Note that this method should be called before `solve`,
    /// once the matrices `H` and `A` have been initialised.
    /// @see solve
    auto decompose(const Matrix& invH, const Matrix& A) -> void;

    /// Solve the KKT problem using an efficient rangespace decomposition approach.
    /// Note that this method requires `decompose` to be called a priori.
    /// @param rhs The right-hand side vector of the KKT equation
    /// @param result[out] The result of the KKT equation
    /// @see decompose
    auto solve(const Vector& a, const Vector& b, Vector& x, Vector& y) -> void;
};

struct KktProblemDiagonalH
{
    /// The internal data for the KKT problem
    Vector invH;
    Matrix AinvH;
    Matrix AinvHAt;
    LLT<Matrix> llt_AinvHAt;

    /// Decompose any necessary matrix before the KKT calculation.
    /// Note that this method should be called before `solve`,
    /// once the matrices `H` and `A` have been initialised.
    /// @see solve
    auto decompose(const Vector& H, const Matrix& A) -> void;

    /// Solve the KKT problem using an efficient rangespace decomposition approach.
    /// Note that this method requires `decompose` to be called a priori.
    /// @param rhs The right-hand side vector of the KKT equation
    /// @param result[out] The result of the KKT equation
    /// @see decompose
    auto solve(const Vector& a, const Vector& b, Vector& x, Vector& y) -> void;
};

struct KktProblemConstantA
{
    /// The pointer to the `H` matrix
    const Matrix* H_ptr;

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
    /// once the matrix `A` has been initialised.
    /// @see solve
    auto decompose(const Matrix& H) -> void;

    /// Solve the KKT problem using an efficient nullspace decomposition approach.
    /// Note that this method requires `decompose` to be called a priori.
    /// @param rhs The right-hand side vector of the KKT equation
    /// @param result[out] The result of the KKT equation
    /// @see decompose
    auto solve(const Vector& a, const Vector& b, Vector& x, Vector& y) -> void;
};

auto KktProblemFull::decompose(const Matrix& H, const Matrix& A) -> void
{
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

auto KktProblemFull::solve(const Vector& a, const Vector& b, Vector& x, Vector& y) -> void
{
    // The dimensions of the KKT problem
    const unsigned n = a.rows();
    const unsigned m = b.rows();

    // Assemble the right-hand side of the KKT equation
    kkt_rhs.segment(0, n) = a;
    kkt_rhs.segment(n, m) = b;

    // Check if the LU decomposition has already been performed
    if(kkt_lu.rows() != n + m or kkt_lu.cols() != n + m)
        error("Cannot solve the KKT equation using the fullspace algorithm.",
            "The LU decomposition of the KKT matrix has not been initialised"
            "or updated for a new problem with different dimension.");

    // Solve the linear system with the LU decomposition already calculated
    kkt_sol = kkt_lu.solve(kkt_rhs);

    // Extract the solution `x` and `y` from the linear system solution `sol`
    x = rows(kkt_sol, 0, n);
    y = rows(kkt_sol, n, m);
}

auto KktProblemInverseH::decompose(const Matrix& invH, const Matrix& A) -> void
{
    invH_ptr = &invH;
    AinvH = A * invH;
    AinvHAt = AinvH * A.transpose();

    llt_AinvHAt.compute(AinvHAt);

    if(llt_AinvHAt.info() != Eigen::Success)
        error("Cannot solve the KKT problem with the rangespace algorithm.",
            "The provided matrix `H` of the KKT equation might not be symmetric positive-definite.");
}

auto KktProblemInverseH::solve(const Vector& a, const Vector& b, Vector& x, Vector& y) -> void
{
    y = llt_AinvHAt.solve(b - AinvH*a);
    x = *invH_ptr * a + AinvH.transpose()*y;
}

auto KktProblemDiagonalH::decompose(const Vector& H, const Matrix& A) -> void
{
    invH = inv(H);
    AinvH = A * diag(invH);
    AinvHAt = AinvH * A.transpose();

    llt_AinvHAt.compute(AinvHAt);

    if(llt_AinvHAt.info() != Eigen::Success)
        error("Cannot solve the KKT problem with the rangespace algorithm.",
            "The provided matrix `H` of the KKT equation might not be symmetric positive-definite.");
}

auto KktProblemDiagonalH::solve(const Vector& a, const Vector& b, Vector& x, Vector& y) -> void
{
    y = llt_AinvHAt.solve(b - AinvH*a);
    x = invH % a + AinvH.transpose()*y;
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

auto KktProblemConstantA::decompose(const Matrix& H) -> void
{
    // Check if the nullspace and rangespace matrices have been initialised
    if(not Z.size())
        error("Cannot solve the KKT equation using the nullspace algorithm.",
            "The matrix `A` of the KKT equation was not initialised.");

    // Set the pointer to the `H` matrix
    H_ptr = &H;

    // Compute the reduced nullspace matrix
    ZtHZ = Z.transpose() * H * Z;

    // Compute Cholesky decomposition of the reduced nullspace matrix
    llt_ZtHZ.compute(ZtHZ);

    // Check if the Cholesky decomposition was successful
    if(llt_ZtHZ.info() != Eigen::Success)
        error("Cannot solve the KKT equation using the nullspace algorithm.",
            "The provided H matrix might not be symmetric positive-definite or is ill-conditioned.");
}

auto KktProblemConstantA::solve(const Vector& a, const Vector& b, Vector& x, Vector& y) -> void
{
    // The dimensions of `x` and `y`
    const unsigned n = x.rows();
    const unsigned m = y.rows();

    // Check if the Cholesky decomposition has already been performed
    if(llt_ZtHZ.rows() != n - m or llt_ZtHZ.cols() != n - m)
        error("Cannot solve the KKT equation using the nullspace algorithm.",
            "The Cholesky decomposition of the reduced Hessian matrix has not "
            "been initialised or updated for a new problem with different dimension.");

    // Compute the `xZ` component of `x`
    xZ = Z.transpose() * (a - *H_ptr*Y*b);
    xZ = llt_ZtHZ.solve(xZ);

    // Compute both `x` and `y` variables
    x = Z*xZ + Y*b;
    y = Y.transpose() * (*H_ptr*x - a);
}

struct KktSolver::Impl
{
    KktScheme scheme = KktScheme::Uninitialised;
    KktProblemFull kkt_dense;
    KktProblemConstantA kkt_const_A;
    KktProblemDiagonalH kkt_diag_H;
    KktProblemInverseH kkt_inv_H;
    KktStatistics statistics;
};

KktSolver::KktSolver()
: pimpl(new Impl())
{}

KktSolver::KktSolver(const KktSolver& other)
: pimpl(new Impl(*other.pimpl))
{}

KktSolver::~KktSolver()
{}

auto KktSolver::operator=(KktSolver other) -> KktSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto KktSolver::statistics() const -> const KktStatistics&
{
    return pimpl->statistics;
}

auto KktSolver::scheme() const -> KktScheme
{
    return pimpl->scheme;
}

auto KktSolver::setConstantA(const Matrix& A) -> void
{
    pimpl->kkt_const_A.initialise(A);
    pimpl->scheme = KktScheme::ConstantA;
}

auto KktSolver::decompose(const Matrix& H, const Matrix& A) -> void
{
    pimpl->kkt_dense.decompose(H, A);
    pimpl->scheme = KktScheme::Full;
}

auto KktSolver::decomposeWithInverseH(const Matrix& invH, const Matrix& A) -> void
{
    pimpl->kkt_inv_H.decompose(invH, A);
    pimpl->scheme = KktScheme::InverseH;
}

auto KktSolver::decomposeWithDiagonalH(const Vector& H, const Matrix& A) -> void
{
    pimpl->kkt_diag_H.decompose(H, A);
    pimpl->scheme = KktScheme::DiagonalH;
}

auto KktSolver::decomposeWithConstantA(const Matrix& H) -> void
{
    if(scheme() != KktScheme::ConstantA)
        error("Cannot perform a KKT decomposition with KktSolver::decomposeWithConstantA.",
              "You have forgotten to call KktSolver::setConstantA once first.");

    pimpl->kkt_const_A.decompose(H);
    pimpl->scheme = KktScheme::ConstantA;
}

auto KktSolver::solve(const Vector& a, const Vector& b, Vector& x, Vector& y) -> void
{
    Time begin = time();

    pimpl->statistics = {};

    switch(scheme())
    {
    case KktScheme::Full: pimpl->kkt_dense.solve(a, b, x, y); break;
    case KktScheme::InverseH: pimpl->kkt_inv_H.solve(a, b, x, y); break;
    case KktScheme::DiagonalH: pimpl->kkt_diag_H.solve(a, b, x, y); break;
    case KktScheme::ConstantA: pimpl->kkt_diag_H.solve(a, b, x, y); break;
    default: error("Cannot solve the KKT equation.",
                   "You have forgotten to call some KktSolver::decompose* "
                   "method before calling KktSolver::solve.");
    }

    pimpl->statistics.succeeded = true;
    pimpl->statistics.time = elapsed(begin);
}

} // namespace Reaktor
