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

// C++ includes
#include <iostream>
// todo remove

// Eigen includes
#include <Reaktor/eigen/Cholesky>
#include <Reaktor/eigen/Core>
#include <Reaktor/eigen/LU>
using namespace Eigen;

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Common/SetUtils.hpp>
#include <Reaktor/Common/TimeUtils.hpp>
#include <Reaktor/Math/MathUtils.hpp>

namespace Reaktor {

struct KktSolverBase
{
    virtual auto decompose(const KktMatrix& lhs) -> void = 0;

    virtual auto solve(const KktVector& rhs, KktSolution& sol) -> void = 0;
};

template<typename LUSolver>
struct KktSolverDense : KktSolverBase
{
    /// The pointer to the left-hand side KKT matrix
    const KktMatrix* lhs;

    /// The internal data for the KKT problem
    Matrix kkt_lhs;
    Vector kkt_rhs;
    Vector kkt_sol;
    LUSolver kkt_lu;

    /// Decompose any necessary matrix before the KKT calculation.
    /// Note that this method should be called before `solve`,
    /// once the matrices `H` and `A` have been initialised.
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
    /// once the matrices `H` and `A` have been initialised.
    virtual auto decompose(const KktMatrix& lhs) -> void;

    /// Solve the KKT problem using an efficient rangespace decomposition approach.
    /// Note that this method requires `decompose` to be called a priori.
    virtual auto solve(const KktVector& rhs, KktSolution& sol) -> void;
};

struct KktSolverRangespaceDiagonal : KktSolverBase
{
    /// The pointer to the left-hand side KKT matrix
    const KktMatrix* lhs;

    /// The matrix `inv(G)` where `G = H + inv(X)*Z`
    Vector invG;
    Matrix AinvG;
    Matrix AinvGAt;
    LLT<Matrix> llt_AinvGAt;

    Vector D;
    Vector H_bar;
    Matrix A_bar;
    Vector a_bar;
    Vector m;

    /// Decompose any necessary matrix before the KKT calculation.
    /// Note that this method should be called before `solve`,
    /// once the matrices `H` and `A` have been initialised.
    virtual auto decompose(const KktMatrix& lhs) -> void;

    /// Solve the KKT problem using an efficient rangespace decomposition approach.
    /// Note that this method requires `decompose` to be called a priori.
    virtual auto solve(const KktVector& rhs, KktSolution& sol) -> void;
};

struct KktSolverRangespaceSparseDiagonal : KktSolverBase
{
    /// The pointer to the left-hand side KKT matrix
    const KktMatrix* lhs;

    Vector D1;
    Matrix A1;
    Matrix A2;
    Vector x1;
    Vector x2;
    Vector z1;
    Vector z2;

    Vector a1, a2;
    Vector c1, c2;

    Vector dx1, dx2, dz1, dz2;

    Indices inonzeros, izeros;

    Vector invH1, invH2;

    Matrix M;
    Matrix m;

    LLT<Matrix> llt_M;

    /// Decompose any necessary matrix before the KKT calculation.
    /// Note that this method should be called before `solve`,
    /// once the matrices `H` and `A` have been initialised.
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

    /// Initialise the constant bottom-left matrix `A` of the KKT equation.
    /// This method should be called once to initialise the `A` matrix
    /// of the KKT equation and whenever it is changed subsequently.
    auto initialise(const Matrix& A) -> void;

    /// Decompose any necessary matrix before the KKT calculation.
    /// Note that this method should be called before `solve`,
    /// once the matrix `A` has been initialised.
    virtual auto decompose(const KktMatrix& lhs) -> void;

    /// Solve the KKT problem using an efficient nullspace decomposition approach.
    /// Note that this method requires `decompose` to be called a priori.
    virtual auto solve(const KktVector& rhs, KktSolution& sol) -> void;
};

template<typename LUSolver>
auto KktSolverDense<LUSolver>::decompose(const KktMatrix& lhs) -> void
{
    /// Update the pointer to the KKT matrix
    this->lhs = &lhs;

    // Check if the Hessian matrix is in the dense mode
    Assert(lhs.H.mode == Hessian::Dense,
        "Cannot solve the KKT equation using PartialPivLU or FullPivLU algorithms.",
        "The Hessian matrix must be in Dense mode.");

    // Auxiliary references to the KKT matrix components
    const auto& x = lhs.x;
    const auto& z = lhs.z;
    const auto& H = lhs.H.dense;
    const auto& A = lhs.A.Ae;

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
auto KktSolverDense<LUSolver>::solve(const KktVector& rhs, KktSolution& sol) -> void
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

    // The dimensions of the KKT problem
    const unsigned n = rx.rows();
    const unsigned m = ry.rows();

    // Assemble the right-hand side of the KKT equation
    kkt_rhs.segment(0, n) = rx + rz/x;
    kkt_rhs.segment(n, m) = ry;

    // Check if the LU decomposition has already been performed
    Assert(kkt_lu.rows() == n + m and kkt_lu.cols() == n + m,
        "Cannot solve the KKT equation using a LU algorithm.",
        "The LU decomposition of the KKT matrix was not performed a priori"
        "or not updated for a new problem with different dimension.");

    // Solve the linear system with the LU decomposition already calculated
    kkt_sol = kkt_lu.solve(kkt_rhs);

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
    const auto& A    = lhs.A.Ae;

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

//auto KktSolverRangespaceDiagonal::decompose(const KktMatrix& lhs) -> void
//{
//    /// Update the pointer to the KKT matrix
//    this->lhs = &lhs;
//
//    // Check if the Hessian matrix is diagonal
//    Assert(lhs.H.mode == Hessian::Diagonal,
//        "Cannot solve the KKT equation using the rangespace algorithm.",
//        "The Hessian matrix must be in Diagonal mode.");
//
//    // Auxiliary references to the KKT matrix components
//    const auto& x = lhs.x;
//    const auto& z = lhs.z;
//    const auto& H = lhs.H.diagonal;
//    const auto& A = lhs.A.Ae;
//
//    invG.noalias() = inv(H + z/x);
//    AinvG.noalias() = A * diag(invG);
//    AinvGAt.noalias() = AinvG * tr(A);
//
//    llt_AinvGAt.compute(AinvGAt);
//
//    Assert(llt_AinvGAt.info() == Eigen::Success,
//        "Cannot solve the KKT problem with the rangespace algorithm.",
//        "The provided Hessian matrix of the KKT equation might not be "
//        "symmetric positive-definite.");
//}
//
//auto KktSolverRangespaceDiagonal::solve(const KktVector& rhs, KktSolution& sol) -> void
//{
//    // Auxiliary references
//    const auto& rx = rhs.rx;
//    const auto& ry = rhs.ry;
//    const auto& rz = rhs.rz;
//    const auto& x  = lhs->x;
//    const auto& z  = lhs->z;
//    auto& dx = sol.dx;
//    auto& dy = sol.dy;
//    auto& dz = sol.dz;
//
//    dy = llt_AinvGAt.solve(ry - AinvG*(rx + rz/x));
//    dx = invG % (rx + rz/x) + tr(AinvG)*dy;
//    dz = (rz - z % dx)/x;
//}

auto KktSolverRangespaceDiagonal::decompose(const KktMatrix& lhs) -> void
{
    /// Update the pointer to the KKT matrix
    this->lhs = &lhs;

    // Check if the Hessian matrix is diagonal
    Assert(lhs.H.mode == Hessian::Diagonal,
        "Cannot solve the KKT equation using the rangespace algorithm.",
        "The Hessian matrix must be in Diagonal mode.");

    // Auxiliary references to the KKT matrix components
    const auto& x = lhs.x;
    const auto& z = lhs.z;
    const auto& H = lhs.H.diagonal;
    const auto& A = lhs.A.Ae;

    D = x.cwiseSqrt();

    H_bar = H % x;
    A_bar = A * diag(D);

    invG.noalias() = inv(H_bar + z);
    AinvG.noalias() = A_bar * diag(invG);
    AinvGAt.noalias() = AinvG * tr(A_bar);

    llt_AinvGAt.compute(AinvGAt);

    Assert(llt_AinvGAt.info() == Eigen::Success,
        "Cannot solve the KKT problem with the rangespace algorithm.",
        "The provided Hessian matrix of the KKT equation might not be "
        "symmetric positive-definite.");
}

auto KktSolverRangespaceDiagonal::solve(const KktVector& rhs, KktSolution& sol) -> void
{
    // Auxiliary references
    const auto& a = rhs.rx;
    const auto& b = rhs.ry;
    const auto& c = rhs.rz;
    const auto& x  = lhs->x;
    const auto& z  = lhs->z;
    auto& dx = sol.dx;
    auto& dy = sol.dy;
    auto& dz = sol.dz;

    a_bar = a%D + c/D;
    m = b - AinvG*a_bar;

    dy = llt_AinvGAt.solve(m);
    dx = invG % a_bar + tr(AinvG)*dy;
    dx = D % dx;
    dz = (c - z % dx)/x;
}

auto KktSolverRangespaceSparseDiagonal::decompose(const KktMatrix& lhs) -> void
{
    /// Update the pointer to the KKT matrix
    this->lhs = &lhs;

    // Check if the Hessian matrix is sparse diagonal
    Assert(lhs.H.mode == Hessian::SparseDiagonal,
        "Cannot decompose the KKT matrix using a sparse diagonal Hessian matrix.",
        "The Hessian matrix was not set to SparseDiagonal mode.");

    // Auxiliary references to the KKT matrix components
    const auto& n = lhs.x.size();
    const auto& x = lhs.x;
    const auto& z = lhs.z;
    const auto& H = lhs.H.sparsediagonal;
    Matrix A = lhs.A.Ae;

    inonzeros = H.inonzeros;
    izeros = range<Index>(n);
    izeros = difference(izeros, inonzeros);

    D1 = H.nonzeros;
    A1 = cols(A, inonzeros);
    A2 = cols(A, izeros);
    x1 = rows(x, inonzeros);
    x2 = rows(x, izeros);
    z1 = rows(z, inonzeros);
    z2 = rows(z, izeros);

    invH1 = inv(D1 + z1/x1);
    invH2 = x2/z2;

    M = A1*diag(invH1)*tr(A1) + A2*diag(invH2)*tr(A2);

    llt_M.compute(M);

    Assert(llt_M.info() == Eigen::Success,
        "Cannot solve the KKT problem with a sparse diagonal Hessian matrix.",
        "The provided sparse diagonal Hessian matrix might not be positive-definite.");
}

auto KktSolverRangespaceSparseDiagonal::solve(const KktVector& rhs, KktSolution& sol) -> void
{
    // Auxiliary references
    const auto& a = rhs.rx;
    const auto& b = rhs.ry;
    const auto& c = rhs.rz;
    auto& dx = sol.dx;
    auto& dy = sol.dy;
    auto& dz = sol.dz;

    a1 = rows(a, inonzeros);
    a2 = rows(a, izeros);
    c1 = rows(c, inonzeros);
    c2 = rows(c, izeros);

    m = b - A1*diag(invH1)*(a1 + c1/x1) - A2*diag(invH2)*(a2 + c2/x2);

//    dy = llt_M.solve(m);
    dy = M.lu().solve(m);
//    dy = M.fullPivLu().solve(m);



    std::cout << "M = \n" << M << std::endl;
    std::cout << "m = \n" << m << std::endl;
    std::cout << "Mdy - m = \n" << (M*dy - m)/m << std::endl;
    std::cout << "dy = \n" << dy << std::endl;

    dz2 = -(a2 + tr(A2)*dy);
    dx2 = (c2 - x2%dz2)/z2;
    dx1 = invH1 % (a1 + c1/x1 + tr(A1)*dy);
//    dx2 = invH2 % (a2 + c2/x2 + tr(A2)*dy);
    dz1 = (c1 - z1%dx1)/x1;

    dx.resize(a.rows());
    dz.resize(c.rows());

    rows(dx, inonzeros) = dx1;
    rows(dz, inonzeros) = dz1;
    rows(dx, izeros)    = dx2;
    rows(dz, izeros)    = dz2;

    Vector res1 = D1%dx1 - tr(A1)*dy - dz1 - a1;
    Vector res2 =        - tr(A2)*dy - dz2 - a2;
    Vector res3 = A1*dx1 + A2*dx2 - b;
    Vector res4 = z1%dx1 + x1%dz1 - c1;
    Vector res5 = z2%dx2 + x2%dz2 - c2;

    std::cout << "res1 = \n" << res1 << std::endl << std::endl;
    std::cout << "res2 = \n" << res2 << std::endl;
    std::cout << "res3 = \n" << res3 << std::endl;
    std::cout << "res4 = \n" << res4 << std::endl;
    std::cout << "res5 = \n" << res5 << std::endl;

//    exit(1);
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

auto KktSolverNullspace::decompose(const KktMatrix& lhs) -> void
{
    /// Update the pointer to the KKT matrix
    this->lhs = &lhs;

    // Check if the Hessian matrix is dense
    Assert(lhs.H.mode == Hessian::Dense,
        "Cannot solve the KKT equation using the nullspace algorithm.",
        "The Hessian matrix must be in the Dense mode.");

    // Auxiliary references to the KKT matrix components
    const auto& x = lhs.x;
    const auto& z = lhs.z;
    const auto& H = lhs.H.dense;
    const auto& A = lhs.A.Ae;

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
    Assert(llt_ZtGZ.rows() == n - m or llt_ZtGZ.cols() == n - m,
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
    KktSolverRangespaceSparseDiagonal kkt_rangespace_sparsediagonal;
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

        if(lhs.H.mode == Hessian::SparseDiagonal)
            base = &kkt_rangespace_sparsediagonal;

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

        if(lhs.H.mode == Hessian::SparseDiagonal)
            base = &kkt_rangespace_sparsediagonal;

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

    result.succeeded = true;
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

} // namespace Reaktor
