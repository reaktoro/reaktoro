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
#include <Reaktor/Math/MathUtils.hpp>
#include <Reaktor/Optimization/OptimumState.hpp>
#include <Reaktor/Optimization/OptimumResult.hpp>

namespace Reaktor {

struct KktSolverBase
{
    virtual auto decompose(const OptimumState& state) -> void = 0;

    virtual auto solve(const Vector& a, const Vector& b, Vector& x, Vector& y) -> void = 0;
};

template<typename LUSolver>
struct KktSolverDense : KktSolverBase
{
    /// The internal data for the KKT problem
    Matrix kkt_lhs;
    Vector kkt_rhs;
    Vector kkt_sol;
    LUSolver kkt_lu;

    /// Decompose any necessary matrix before the KKT calculation.
    /// Note that this method should be called before `solve`,
    /// once the matrices `H` and `A` have been initialised.
    virtual auto decompose(const OptimumState& state) -> void;

    /// Solve the KKT problem using a dense LU decomposition.
    /// Note that this method requires `decompose` to be called a priori.
    virtual auto solve(const Vector& a, const Vector& b, Vector& x, Vector& y) -> void;
};

struct KktSolverRangespaceInverse : KktSolverBase
{
    /// The matrix `inv(G)` where `G = H + inv(X)*Z`
    Matrix invG;
    Matrix AinvG;
    Matrix AinvGAt;
    LLT<Matrix> llt_AinvGAt;

    /// Decompose any necessary matrix before the KKT calculation.
    /// Note that this method should be called before `solve`,
    /// once the matrices `H` and `A` have been initialised.
    virtual auto decompose(const OptimumState& state) -> void;

    /// Solve the KKT problem using an efficient rangespace decomposition approach.
    /// Note that this method requires `decompose` to be called a priori.
    virtual auto solve(const Vector& a, const Vector& b, Vector& x, Vector& y) -> void;
};

struct KktSolverRangespaceDiagonal : KktSolverBase
{
    /// The matrix `inv(G)` where `G = H + inv(X)*Z`
    Vector invG;
    Matrix AinvG;
    Matrix AinvGAt;
    LLT<Matrix> llt_AinvGAt;

    /// Decompose any necessary matrix before the KKT calculation.
    /// Note that this method should be called before `solve`,
    /// once the matrices `H` and `A` have been initialised.
    virtual auto decompose(const OptimumState& state) -> void;

    /// Solve the KKT problem using an efficient rangespace decomposition approach.
    /// Note that this method requires `decompose` to be called a priori.
    virtual auto solve(const Vector& a, const Vector& b, Vector& x, Vector& y) -> void;
};

struct KktSolverNullspace : KktSolverBase
{
    /// The matrix `A` of the KKT problem
    Matrix A;

    /// The matrix `G = H + inv(X)*Z` of the KKT equation
    Matrix G;

    /// The nullspace matrix of `A` with the property `AZ = 0`
    Matrix Z;

    /// The rangespace matrix of `A` with the property `AY = I`
    Matrix Y;

    /// Auxiliary data for the nullspace algorithm
    Matrix ZtGZ;
    LLT<Matrix> llt_ZtGZ;
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
    virtual auto decompose(const OptimumState& state) -> void;

    /// Solve the KKT problem using an efficient nullspace decomposition approach.
    /// Note that this method requires `decompose` to be called a priori.
    virtual auto solve(const Vector& a, const Vector& b, Vector& x, Vector& y) -> void;
};

template<typename LUSolver>
auto KktSolverDense<LUSolver>::decompose(const OptimumState& state) -> void
{
    // Check if the Hessian matrix is in the dense mode
    Assert(state.H.mode == Hessian::Dense,
        "Cannot solve the KKT equation using PartialPivLU or FullPivLU algorithms.",
        "The Hessian matrix must be in Dense mode.");

    // Define some auxiliary references to variables
    const auto& x = state.x;
    const auto& z = state.z;
    const auto& H = state.H.dense;
    const auto& A = state.A;

    // The dimensions of the KKT problem
    const unsigned n = A.cols();
    const unsigned m = A.rows();

    // Ensure the components of the KKT equation have adequate dimensions
    kkt_lhs.resize(n + m, n + m);
    kkt_rhs.resize(n + m);
    kkt_sol.resize(n + m);

    // Assemble the left-hand side of the KKT equation
    kkt_lhs.block(0, 0, n, n).noalias()   =  H;
    kkt_lhs.block(0, 0, n, n).diagonal() +=  z/x;
    kkt_lhs.block(0, n, n, m).noalias()   = -tr(A);
    kkt_lhs.block(n, 0, m, n).noalias()   =  A;
    kkt_lhs.block(n, n, m, m).noalias()   =  zeros(m, m);

    // Perform the LU decomposition
    kkt_lu.compute(kkt_lhs);
}

template<typename LUSolver>
auto KktSolverDense<LUSolver>::solve(const Vector& a, const Vector& b, Vector& x, Vector& y) -> void
{
    // The dimensions of the KKT problem
    const unsigned n = a.rows();
    const unsigned m = b.rows();

    // Assemble the right-hand side of the KKT equation
    kkt_rhs.segment(0, n) = a;
    kkt_rhs.segment(n, m) = b;

    // Check if the LU decomposition has already been performed
    Assert(kkt_lu.rows() == n + m and kkt_lu.cols() == n + m,
        "Cannot solve the KKT equation using a LU algorithm.",
        "The LU decomposition of the KKT matrix was not performed a priori"
        "or not updated for a new problem with different dimension.");

    // Solve the linear system with the LU decomposition already calculated
    kkt_sol = kkt_lu.solve(kkt_rhs);

    // Extract the solution `x` and `y` from the linear system solution `sol`
    x = rows(kkt_sol, 0, n);
    y = rows(kkt_sol, n, m);
}

auto KktSolverRangespaceInverse::decompose(const OptimumState& state) -> void
{
    // Check if the Hessian matrix is in inverse more
    Assert(state.H.mode == Hessian::Inverse,
        "Cannot solve the KKT equation using the rangespace algorithm.",
        "The Hessian matrix must be in Diagonal or Inverse mode.");

    // Define some auxiliary references to variables
    const auto& x = state.x;
    const auto& z = state.z;
    const auto& invH = state.H.inverse;
    const auto& A = state.A;

    invG = inverseShermanMorrison(invH, z/x);
    AinvG = A * invG;
    AinvGAt = AinvG * tr(A);

    llt_AinvGAt.compute(AinvGAt);

    Assert(llt_AinvGAt.info() == Eigen::Success,
        "Cannot solve the KKT problem with the rangespace algorithm.",
        "The provided Hessian matrix of the KKT equation might not be "
        "symmetric positive-definite.");
}

auto KktSolverRangespaceInverse::solve(const Vector& a, const Vector& b, Vector& x, Vector& y) -> void
{
    y = llt_AinvGAt.solve(b - AinvG*a);
    x = invG * a + tr(AinvG)*y;
}

auto KktSolverRangespaceDiagonal::decompose(const OptimumState& state) -> void
{
    // Check if the Hessian matrix is diagonal
    Assert(state.H.mode == Hessian::Diagonal,
        "Cannot solve the KKT equation using the rangespace algorithm.",
        "The Hessian matrix must be in Diagonal or Inverse mode.");

    // Define some auxiliary references to variables
    const auto& x = state.x;
    const auto& z = state.z;
    const auto& H = state.H.diagonal;
    const auto& A = state.A;

    invG.noalias() = inv(H + z/x);
    AinvG.noalias() = A * diag(invG);
    AinvGAt.noalias() = AinvG * tr(A);

    llt_AinvGAt.compute(AinvGAt);

    Assert(llt_AinvGAt.info() == Eigen::Success,
        "Cannot solve the KKT problem with the rangespace algorithm.",
        "The provided Hessian matrix of the KKT equation might not be "
        "symmetric positive-definite.");
}

auto KktSolverRangespaceDiagonal::solve(const Vector& a, const Vector& b, Vector& x, Vector& y) -> void
{
    y = llt_AinvGAt.solve(b - AinvG*a);
    x = invG % a + tr(AinvG)*y;
}

auto KktSolverNullspace::initialise(const Matrix& newA) -> void
{
    // Check if `newA` was used last time to avoid repeated operations
    if(A.rows() == newA.rows() and
       A.cols() == newA.cols() and
       A == newA) return;

    // Update the matrix `A` of the KKT equation
    A.noalias() = newA;

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

auto KktSolverNullspace::decompose(const OptimumState& state) -> void
{
    // Check if the Hessian matrix is dense
    Assert(state.H.mode == Hessian::Dense,
        "Cannot solve the KKT equation using the nullspace algorithm.",
        "The Hessian matrix must be in the Dense mode.");

    // Define some auxiliary references to variables
    const auto& x = state.x;
    const auto& z = state.z;
    const auto& H = state.H.dense;
    const auto& A = state.A;

    // Initialise the solver with the matrix `A`
    initialise(A);

    // Set matrix `G = H + inv(X)*Z`
    G.noalias() = H;
    G.diagonal() += z/x;

    // Compute the reduced Hessian matrix
    ZtGZ = tr(Z) * G * Z;

    // Compute Cholesky decomposition of the reduced Hessian matrix
    llt_ZtGZ.compute(ZtGZ);

    // Check if the Cholesky decomposition was successful
    Assert(llt_ZtGZ.info() == Eigen::Success,
        "Cannot solve the KKT equation using the nullspace algorithm.",
        "The provided H matrix might not be symmetric positive-definite or "
        "is ill-conditioned.");
}

auto KktSolverNullspace::solve(const Vector& a, const Vector& b, Vector& x, Vector& y) -> void
{
    // The dimensions of `x` and `y`
    const unsigned n = x.rows();
    const unsigned m = y.rows();

    // Check if the Cholesky decomposition has already been performed
    Assert(llt_ZtGZ.rows() == n - m or llt_ZtGZ.cols() == n - m,
        "Cannot solve the KKT equation using the nullspace algorithm.",
        "The Cholesky decomposition of the reduced Hessian matrix has not "
        "been performed a priori or not updated for a new problem with "
        "different dimension.");

    // Compute the `xZ` component of `x`
    xZ = Z.transpose() * (a - G*Y*b);
    xZ = llt_ZtGZ.solve(xZ);

    // Compute both `x` and `y` variables
    x = Z*xZ + Y*b;
    y = Y.transpose() * (G*x - a);
}

struct KktSolver::Impl
{
    KktInfo info;
    KktOptions options;
    KktSolverDense<PartialPivLU<Matrix>> kkt_partial_lu;
    KktSolverDense<FullPivLU<Matrix>> kkt_full_lu;
    KktSolverNullspace kkt_nullspace;
    KktSolverRangespaceDiagonal kkt_rangespace_diagonal;
    KktSolverRangespaceInverse kkt_rangespace_inverse;
    KktSolverBase* base;

    auto decompose(const OptimumState& state) -> void;

    auto solve(const Vector& a, const Vector& b, Vector& x, Vector& y) -> void;
};

auto KktSolver::Impl::decompose(const OptimumState& state) -> void
{
    if(options.method == KktMethod::Automatic)
    {
        if(state.H.mode == Hessian::Dense)
            base = &kkt_partial_lu;

        if(state.H.mode == Hessian::Diagonal)
            base = &kkt_rangespace_diagonal;

        if(state.H.mode == Hessian::Inverse)
            base = &kkt_rangespace_inverse;
    }

    if(options.method == KktMethod::FullPivLU)
        base = &kkt_full_lu;

    if(options.method == KktMethod::Nullspace)
        base = &kkt_nullspace;

    if(options.method == KktMethod::Rangespace)
    {
        if(state.H.mode == Hessian::Diagonal)
            base = &kkt_rangespace_diagonal;

        if(state.H.mode == Hessian::Inverse)
            base = &kkt_rangespace_inverse;
    }

    Time begin = time();

    base->decompose(state);

    info.succeeded = true;
    info.time_decompose = elapsed(begin);
}

auto KktSolver::Impl::solve(const Vector& a, const Vector& b, Vector& x, Vector& y) -> void
{
    Time begin = time();

    base->solve(a, b, x, y);

    info.succeeded = true;
    info.time_solve = elapsed(begin);
}

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

auto KktSolver::info() const -> const KktInfo&
{
    return pimpl->info;
}

auto KktSolver::setOptions(const KktOptions& options) -> void
{
    pimpl->options = options;
}

auto KktSolver::decompose(const OptimumState& state) -> void
{
    pimpl->decompose(state);
}

auto KktSolver::solve(const Vector& a, const Vector& b, Vector& dx, Vector& dy) -> void
{
    pimpl->solve(a, b, dx, dy);
}

} // namespace Reaktor
