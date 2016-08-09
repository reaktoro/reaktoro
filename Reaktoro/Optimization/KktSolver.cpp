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

#include "KktSolver.hpp"

// Eigen includes
#include <Reaktoro/Eigen/LU>
#include <Reaktoro/Eigen/Cholesky>
using namespace Eigen;

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/TimeUtils.hpp>
#include <Reaktoro/Math/MathUtils.hpp>

namespace Reaktoro {

struct KktSolverBase
{
    virtual auto decompose(const KktMatrix& lhs) -> void = 0;

    virtual auto solve(const KktVector& rhs, KktSolution& sol) -> void = 0;
};

template<typename LUSolver>
struct KktSolverDense : KktSolverBase
{
    /// The vectors x and z
    Vector x, z;

    /// The internal data for the KKT problem
    Matrix kkt_lhs;
    Vector kkt_rhs;
    Vector kkt_sol;
    LUSolver kkt_lu;

    /// Decompose any necessary matrix before the KKT calculation.
    /// Note that this method should be called before `solve`,
    /// once the matrices `H` and `A` have been initialized.
    virtual auto decompose(const KktMatrix& lhs) -> void;

    /// Solve the KKT problem using a dense LU decomposition.
    /// Note that this method requires `decompose` to be called a priori.
    virtual auto solve(const KktVector& rhs, KktSolution& sol) -> void;
};

struct KktSolverRangespaceInverse : KktSolverBase
{
    /// The pointer to the left-hand side KKT matrix
    const KktMatrix* lhs;

    /// The matrix `inv(G)` where `G = H + inv(X)*Z`
    Matrix invG;
    Matrix AinvG;
    Matrix AinvGAt;
    LLT<Matrix> llt_AinvGAt;

    /// Decompose any necessary matrix before the KKT calculation.
    /// Note that this method should be called before `solve`,
    /// once the matrices `H` and `A` have been initialized.
    virtual auto decompose(const KktMatrix& lhs) -> void;

    /// Solve the KKT problem using an efficient rangespace decomposition approach.
    /// Note that this method requires `decompose` to be called a priori.
    virtual auto solve(const KktVector& rhs, KktSolution& sol) -> void;
};

struct KktSolverRangespaceDiagonal : KktSolverBase
{
    Indices ipivot, inonpivot;

    Vector X, Z;
    Vector D, D1, D2;
    Matrix A1, A2;
    Vector a1, a2;
    Vector dx1, dx2;
    Vector r;

    Vector invD1;
    Matrix A1invD1;
    Matrix A1invD1A1t;

    Vector kkt_rhs, kkt_sol;
    Matrix kkt_lhs;

    PartialPivLU<Matrix> lu;

    /// Decompose any necessary matrix before the KKT calculation.
    /// Note that this method should be called before `solve`,
    /// once the matrices `H` and `A` have been initialized.
    virtual auto decompose(const KktMatrix& lhs) -> void;

    /// Solve the KKT problem using an efficient rangespace decomposition approach.
    /// Note that this method requires `decompose` to be called a priori.
    virtual auto solve(const KktVector& rhs, KktSolution& sol) -> void;
};

struct KktSolverNullspace : KktSolverBase
{
    /// The pointer to the left-hand side KKT matrix
    const KktMatrix* lhs;

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

    /// Initialize the constant bottom-left matrix `A` of the KKT equation.
    /// This method should be called once to initialize the `A` matrix
    /// of the KKT equation and whenever it is changed subsequently.
    auto initialize(const Matrix& A) -> void;

    /// Decompose any necessary matrix before the KKT calculation.
    /// Note that this method should be called before `solve`,
    /// once the matrix `A` has been initialized.
    virtual auto decompose(const KktMatrix& lhs) -> void;

    /// Solve the KKT problem using an efficient nullspace decomposition approach.
    /// Note that this method requires `decompose` to be called a priori.
    virtual auto solve(const KktVector& rhs, KktSolution& sol) -> void;
};

template<typename LUSolver>
auto KktSolverDense<LUSolver>::decompose(const KktMatrix& lhs) -> void
{
    /// Update x and z
    x = lhs.x;
    z = lhs.z;

    // Check if the Hessian matrix is in the dense mode
    Assert(lhs.H.mode == Hessian::Dense || lhs.H.mode == Hessian::Diagonal,
        "Cannot solve the KKT equation using PartialPivLU or FullPivLU algorithms.",
        "The Hessian matrix must be in Dense or Diagonal mode.");

    // Auxiliary references to the KKT matrix components
    const auto& H = lhs.H;
    const auto& A = lhs.A;
    const auto& gamma = lhs.gamma;
    const auto& delta = lhs.delta;

    // The dimensions of the KKT problem
    const unsigned n = A.cols();
    const unsigned m = A.rows();

    // Ensure the components of the KKT equation have adequate dimensions
    kkt_lhs.resize(n + m, n + m);
    kkt_rhs.resize(n + m);
    kkt_sol.resize(n + m);

    // Assemble the left-hand side of the KKT equation
    if(H.mode == Hessian::Dense) kkt_lhs.block(0, 0, n, n).noalias() = H.dense;
    else kkt_lhs.block(0, 0, n, n) = diag(H.diagonal);
    kkt_lhs.block(0, 0, n, n).diagonal() +=  z/x;
    kkt_lhs.block(0, 0, n, n).diagonal() +=  gamma*gamma*ones(n);
    kkt_lhs.block(0, n, n, m).noalias() = -tr(A);
    kkt_lhs.block(n, 0, m, n).noalias() =  A;
    kkt_lhs.block(n, n, m, m).noalias() = delta*delta*identity(m, m);

    // Perform the LU decomposition
    kkt_lu.compute(kkt_lhs);
}

template<typename LUSolver>
auto KktSolverDense<LUSolver>::solve(const KktVector& rhs, KktSolution& sol) -> void
{
    // Auxiliary references
    const auto& rx = rhs.rx;
    const auto& ry = rhs.ry;
    const auto& rz = rhs.rz;
    auto& dx = sol.dx;
    auto& dy = sol.dy;
    auto& dz = sol.dz;

    // The dimensions of the KKT problem
    const unsigned n = rx.rows();
    const unsigned m = ry.rows();

    // Assemble the right-hand side of the KKT equation
    kkt_rhs.segment(0, n) = rx + rz/x;
    kkt_rhs.segment(n, m) = ry;

    // Check if the LU decomposition has already been performed
    Assert(kkt_lu.rows() == n + m && kkt_lu.cols() == n + m,
        "Cannot solve the KKT equation using a LU algorithm.",
        "The LU decomposition of the KKT matrix was not performed a priori"
        "or not updated for a new problem with different dimension.");

    // Solve the linear system with the LU decomposition already calculated
    kkt_sol = kkt_lu.solve(kkt_rhs);

    // If the solution failed before (perhaps because PartialPivLU was used), use FullPivLU
    if(!kkt_sol.allFinite())
        kkt_sol = kkt_lhs.fullPivLu().solve(kkt_rhs);

    // Extract the solution `x` and `y` from the linear system solution `sol`
    dx = rows(kkt_sol, 0, n);
    dy = rows(kkt_sol, n, m);
    dz = (rz - z % dx)/x;
}

auto KktSolverRangespaceInverse::decompose(const KktMatrix& lhs) -> void
{
    /// Update the pointer to the KKT matrix
    this->lhs = &lhs;

    // Check if the Hessian matrix is in inverse more
    Assert(lhs.H.mode == Hessian::Inverse,
        "Cannot solve the KKT equation using the rangespace algorithm.",
        "The Hessian matrix must be in Inverse mode.");

    // Auxiliary references to the KKT matrix components
    const auto& x    = lhs.x;
    const auto& z    = lhs.z;
    const auto& invH = lhs.H.inverse;
    const auto& A    = lhs.A;

    invG = inverseShermanMorrison(invH, z/x);
    AinvG = A * invG;
    AinvGAt = AinvG * tr(A);

    llt_AinvGAt.compute(AinvGAt);

    Assert(llt_AinvGAt.info() == Eigen::Success,
        "Cannot solve the KKT problem with the rangespace algorithm.",
        "The provided Hessian matrix of the KKT equation might not be "
        "symmetric positive-definite.");
}

auto KktSolverRangespaceInverse::solve(const KktVector& rhs, KktSolution& sol) -> void
{
    // Auxiliary references
    const auto& rx = rhs.rx;
    const auto& ry = rhs.ry;
    const auto& rz = rhs.rz;
    const auto& x  = lhs->x;
    const auto& z  = lhs->z;
    auto& dx = sol.dx;
    auto& dy = sol.dy;
    auto& dz = sol.dz;

    dy = llt_AinvGAt.solve(ry - AinvG*(rx + rz/x));
    dx = invG * (rx + rz/x) + tr(AinvG)*dy;
    dz = (rz - z % dx)/x;
}

auto KktSolverRangespaceDiagonal::decompose(const KktMatrix& lhs) -> void
{
    // Check if the Hessian matrix is diagonal
    Assert(lhs.H.mode == Hessian::Diagonal,
        "Cannot solve the KKT equation using the rangespace algorithm.",
        "The Hessian matrix must be in Diagonal mode.");

    // Initialize diagonal matrices X and Z
    X = lhs.x;
    Z = lhs.z;

    // Auxiliary references to the KKT matrix components
    const auto& A = lhs.A;
    const auto& H = lhs.H.diagonal;
    const auto& gamma = lhs.gamma;
    const auto& delta = lhs.delta;

    const unsigned n = A.cols();
    const unsigned m = A.rows();

    D.noalias() = H + Z/X + gamma*gamma*ones(n);

    ipivot.clear();
    inonpivot.clear();
    ipivot.reserve(n);
    inonpivot.reserve(n);
    for(unsigned i = 0; i < n; ++i)
        if(D[i] > norminf(A.col(i))) ipivot.push_back(i);
        else inonpivot.push_back(i);

    rows(D, ipivot).to(D1);
    rows(D, inonpivot).to(D2);
    cols(A, ipivot).to(A1);
    cols(A, inonpivot).to(A2);

    invD1.noalias() = inv(D1);
    A1invD1.noalias() = A1*diag(invD1);
    A1invD1A1t.noalias() = A1invD1*tr(A1);

    const unsigned n2 = inonpivot.size();
    const unsigned t  = m + n2;

    kkt_lhs = zeros(t, t);
    kkt_lhs.topLeftCorner(n2, n2).diagonal() = D2;
    kkt_lhs.topRightCorner(n2, m).noalias() = -tr(A2);
    kkt_lhs.bottomLeftCorner(m, n2).noalias() = A2;
    kkt_lhs.bottomRightCorner(m, m).noalias() = A1invD1A1t;
    kkt_lhs.bottomRightCorner(m, m).diagonal() += delta*delta*ones(m);

    lu.compute(kkt_lhs);
}

auto KktSolverRangespaceDiagonal::solve(const KktVector& rhs, KktSolution& sol) -> void
{
    // Auxiliary references
    const auto& a = rhs.rx;
    const auto& b = rhs.ry;
    const auto& c = rhs.rz;
    auto& dx = sol.dx;
    auto& dy = sol.dy;
    auto& dz = sol.dz;

    r.noalias() = a + c/X;
    rows(r, ipivot).to(a1);
    rows(r, inonpivot).to(a2);

    const unsigned n1 = A1.cols();
    const unsigned n2 = A2.cols();
    const unsigned n  = n1 + n2;
    const unsigned m  = A1.rows();
    const unsigned t  = n2 + m;

    kkt_rhs.resize(t);
    kkt_rhs.segment( 0, n2).noalias() = a2;
    kkt_rhs.segment(n2,  m).noalias() = b - A1invD1*a1;

    kkt_sol.noalias() = lu.solve(kkt_rhs);

    if(!kkt_sol.allFinite())
        kkt_sol = kkt_lhs.fullPivLu().solve(kkt_rhs);

    dy.noalias() = kkt_sol.segment(n2, m);

    dx1.noalias() = a1 % invD1 + tr(A1invD1)*dy;
    dx2.noalias() = kkt_sol.segment(0, n2);

    dx.resize(n);
    rows(dx, ipivot)    = dx1;
    rows(dx, inonpivot) = dx2;

    dz.noalias() = (c - Z % dx)/X;
}

auto KktSolverNullspace::initialize(const Matrix& newA) -> void
{
    // Check if `newA` was used last time to avoid repeated operations
    if(A.rows() == newA.rows() &&
       A.cols() == newA.cols() &&
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

auto KktSolverNullspace::decompose(const KktMatrix& lhs) -> void
{
    /// Update the pointer to the KKT matrix
    this->lhs = &lhs;

    // Check if the Hessian matrix is dense
    Assert(lhs.H.mode == Hessian::Dense || lhs.H.mode == Hessian::Diagonal,
        "Cannot solve the KKT equation using the nullspace algorithm.",
        "The Hessian matrix must be either in the Dense or Diagonal mode.");

    // Auxiliary references to the KKT matrix components
    const auto& x = lhs.x;
    const auto& z = lhs.z;
    const auto& H = lhs.H;
    const auto& A = lhs.A;

    // Initialize the solver with the matrix `A`
    initialize(A);

    // Set matrix `G = H + inv(X)*Z`
    G.noalias() = (H.mode == Hessian::Dense) ? H.dense : diag(H.diagonal);
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

auto KktSolverNullspace::solve(const KktVector& rhs, KktSolution& sol) -> void
{
    // Auxiliary references
    const auto& rx = rhs.rx;
    const auto& ry = rhs.ry;
    const auto& rz = rhs.rz;
    const auto& x  = lhs->x;
    const auto& z  = lhs->z;
    auto& dx = sol.dx;
    auto& dy = sol.dy;
    auto& dz = sol.dz;

    // The dimensions of `x` and `y`
    const unsigned n = rx.rows();
    const unsigned m = ry.rows();

    // Check if the Cholesky decomposition has already been performed
    Assert(llt_ZtGZ.rows() == n - m || llt_ZtGZ.cols() == n - m,
        "Cannot solve the KKT equation using the nullspace algorithm.",
        "The Cholesky decomposition of the reduced Hessian matrix has not "
        "been performed a priori or not updated for a new problem with "
        "different dimension.");

    // Compute the `xZ` component of `x`
    xZ = Z.transpose() * ((rx + rz/x) - G*Y*ry);
    xZ = llt_ZtGZ.solve(xZ);

    // Compute both `x` and `y` variables
    dx = Z*xZ + Y*ry;
    dy = Y.transpose() * (G*x - (rx + rz/x));
    dz = (rz - z % dx)/x;
}

struct KktSolver::Impl
{
    KktResult result;
    KktOptions options;
    KktSolverDense<PartialPivLU<Matrix>> kkt_partial_lu;
    KktSolverDense<FullPivLU<Matrix>> kkt_full_lu;
    KktSolverNullspace kkt_nullspace;
    KktSolverRangespaceDiagonal kkt_rangespace_diagonal;
    KktSolverRangespaceInverse kkt_rangespace_inverse;
    KktSolverBase* base;

    auto decompose(const KktMatrix& lhs) -> void;

    auto solve(const KktVector& rhs, KktSolution& sol) -> void;
};

auto KktSolver::Impl::decompose(const KktMatrix& lhs) -> void
{
    if(options.method == KktMethod::Automatic)
    {
        if(lhs.H.mode == Hessian::Dense)
            base = &kkt_partial_lu;

        if(lhs.H.mode == Hessian::Diagonal)
            base = &kkt_rangespace_diagonal;

        if(lhs.H.mode == Hessian::Inverse)
            base = &kkt_rangespace_inverse;
    }

    if(options.method == KktMethod::PartialPivLU)
        base = &kkt_partial_lu;

    if(options.method == KktMethod::FullPivLU)
        base = &kkt_full_lu;

    if(options.method == KktMethod::Nullspace)
        base = &kkt_nullspace;

    if(options.method == KktMethod::Rangespace)
    {
        if(lhs.H.mode == Hessian::Diagonal)
            base = &kkt_rangespace_diagonal;

        if(lhs.H.mode == Hessian::Inverse)
            base = &kkt_rangespace_inverse;
    }

    Time begin = time();

    base->decompose(lhs);

    result.succeeded = true;
    result.time_decompose = elapsed(begin);
}

auto KktSolver::Impl::solve(const KktVector& rhs, KktSolution& sol) -> void
{
    Time begin = time();

    base->solve(rhs, sol);

    result.succeeded = sol.dx.allFinite() && sol.dy.allFinite() && sol.dz.allFinite();
    result.time_solve = elapsed(begin);
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

auto KktSolver::result() const -> const KktResult&
{
    return pimpl->result;
}

auto KktSolver::setOptions(const KktOptions& options) -> void
{
    pimpl->options = options;
}

auto KktSolver::decompose(const KktMatrix& lhs) -> void
{
    pimpl->decompose(lhs);
}

auto KktSolver::solve(const KktVector& rhs, KktSolution& sol) -> void
{
    pimpl->solve(rhs, sol);
}

} // namespace Reaktoro
