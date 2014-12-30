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

struct KktSolver::Impl
{
    Impl();

    auto initialise(const KktProblem& problem, const KktOptions& options) -> void;

    auto computeZY(const KktProblem& problem, const KktOptions& options) -> void;

    auto scaleZY(const KktProblem& problem, const KktOptions& options) -> void;

    auto solveFullspaceDense(const KktProblem& problem, KktResult& result, const KktOptions& options) -> void;

    auto solveFullspaceSparse(const KktProblem& problem, KktResult& result, const KktOptions& options) -> void;

    auto solveRangespace(const KktProblem& problem, KktResult& result, const KktOptions& options) -> void;

    auto solveNullspace(const KktProblem& problem, KktResult& result, const KktOptions& options) -> void;

    auto solveNullspacePartial(const KktProblem& problem, KktResult& result, const KktOptions& options) -> void;

    auto solve(const KktProblem& problem, KktResult& result, const KktOptions& options) -> void;

    bool has_any_scaling = false;
    bool has_xscaling = false;
    bool has_yscaling = false;

    KktProblem scaled_problem;

    /// Nulspace related variables
    bool computed_ZY_before = false;
    Matrix Z;
    Matrix Y;
    Matrix scaled_Z;
    Matrix scaled_Y;
    Matrix ZtHZ;
    Vector xZ;

    /// Rangespace related variables
    Vector invH_diag;
    Matrix AinvH;
    Matrix AinvHAt;
    LLT<Matrix> llt_AinvHAt;

    /// Fullspace related variables
    Matrix lhs;
    Vector rhs;
    Vector sol;
};

KktSolver::Impl::Impl()
{}

auto KktSolver::Impl::initialise(const KktProblem& problem, const KktOptions& options) -> void
{
    // Check if there is scaling on `x` or `y`
    has_xscaling = options.xscaling.size();
    has_yscaling = options.yscaling.size();

    // Check if there is any scaling on `x` or `y`
    has_any_scaling = has_xscaling or has_yscaling;

    // The diagonal representation of the scaling matrices for `x` and `y`
    const auto& dx = options.xscaling;
    const auto& dy = options.yscaling;
    const auto Dx = diag(dx);
    const auto Dy = diag(dy);

    // Initialise the scaled KKT members
    if(has_xscaling and has_yscaling)
    {
        if(options.diagonalH)
            scaled_problem.H = diag(dx % diagonal(problem.H) % dx);
        else
            scaled_problem.H = Dx * problem.H * Dx;
        scaled_problem.A = Dy * problem.A * Dx;
        scaled_problem.f = Dx * problem.f;
        scaled_problem.g = Dy * problem.g;
    }
    else if(has_xscaling)
    {
        if(options.diagonalH)
            scaled_problem.H = diag(dx % diagonal(problem.H) % dx);
        else
            scaled_problem.H = Dx * problem.H * Dx;
        scaled_problem.A = problem.A * Dx;
        scaled_problem.f = Dx * problem.f;
        scaled_problem.g = problem.g;
    }
    else if(has_yscaling)
    {
        scaled_problem.H = problem.H;
        scaled_problem.A = Dy * problem.A;
        scaled_problem.f = problem.f;
        scaled_problem.g = Dy * problem.g;
    }
}

auto KktSolver::Impl::computeZY(const KktProblem& problem, const KktOptions& options) -> void
{
    // Check if the nullspace and rangespace matrices `Z` and `Y` need to be updated
    if(not computed_ZY_before or not options.persistentA)
    {
        // Update the data-member `computed_ZY_before`
        computed_ZY_before = true;

        // The dimensions of the matrix `A`
        const unsigned n = problem.A.cols();
        const unsigned m = problem.A.rows();

        // Perform a LU decomposition of the matrix `A`
        const FullPivLU<Matrix> lu_A = problem.A.fullPivLu();
        const Matrix L = lu_A.matrixLU().leftCols(m).triangularView<UnitLower>();
        const Matrix U = lu_A.matrixLU().triangularView<Upper>();
        const auto P1 = lu_A.permutationP();
        const auto P2 = lu_A.permutationQ();
        const auto U1 = U.leftCols(m);
        const auto U2 = U.rightCols(n - m);

        // Update the data-member nullspace matrix `Z` of `A`
        Z = zeros(n, n - m);
        Z.topRows(m) = -U1.triangularView<Upper>().solve(U2);
        Z.bottomRows(n - m) = identity(n - m, n - m);
        Z = P2*Z;

        // Update the data-member rangespace matrix `Y` of `A`
        Y = zeros(n, m);
        Y.topRows(m) = L.triangularView<Lower>().solve(identity(m, m));
        Y.topRows(m) = U1.triangularView<Upper>().solve(Y.topRows(m));
        Y = P2*Y*P1;
    }

    scaleZY(problem, options);
}

auto KktSolver::Impl::scaleZY(const KktProblem& problem, const KktOptions& options) -> void
{
    // Compute the scaled versions of `Z` and `Y`
    if(has_xscaling and has_yscaling)
    {
        // The diagonal representation of the scaling matrices for `x` and `y`
        const auto inv_Dx = options.xscaling.cwiseInverse().asDiagonal();
        const auto inv_Dy = options.yscaling.cwiseInverse().asDiagonal();
        scaled_Z.noalias() = inv_Dx * Z;
        scaled_Y.noalias() = inv_Dx * Y * inv_Dy;
    }
    else if(has_xscaling)
    {
        // The diagonal representation of the scaling matrix for `x`
        const auto inv_Dx = options.xscaling.cwiseInverse().asDiagonal();
        scaled_Z.noalias() = inv_Dx * Z;
        scaled_Y.noalias() = inv_Dx * Y;
    }
    else if(has_yscaling)
    {
        // The diagonal representation of the scaling matrix for `y`
        const auto inv_Dy = options.yscaling.cwiseInverse().asDiagonal();
        scaled_Z.noalias() = Z;
        scaled_Y.noalias() = Y * inv_Dy;
    }
}

auto KktSolver::Impl::solveFullspaceDense(const KktProblem& problem, KktResult& result, const KktOptions& options) -> void
{
    // The references to the components of the KKT problem
    const auto& A = has_any_scaling ? scaled_problem.A : problem.A;
    const auto& H = has_any_scaling ? scaled_problem.H : problem.H;
    const auto& f = has_any_scaling ? scaled_problem.f : problem.f;
    const auto& g = has_any_scaling ? scaled_problem.g : problem.g;

    // The dimensions of the KKT problem
    const unsigned n = A.cols();
    const unsigned m = A.rows();

    // Ensure the components of the KKT equation have adequate dimensions
    lhs.resize(n + m, n + m);
    rhs.resize(n + m);
    sol.resize(n + m);

    // Assemble the left-hand side of the KKT equation
    lhs.block(0, 0, n, n).noalias() =  H;
    lhs.block(0, n, n, m).noalias() = -A.transpose();
    lhs.block(n, 0, m, n).noalias() =  A;
    lhs.block(n, n, m, m).noalias() =  zeros(m, m);

    // Assemble the right-hand side of the KKT equation
    rhs.segment(0, n) = f;
    rhs.segment(n, m) = g;

    // Solve with full partial pivoting for higher stability
    sol = lhs.fullPivLu().solve(rhs);

    // Extract the solution `x` and `y` from the linear system solution `sol`
    result.solution.x = rows(sol, 0, n);
    result.solution.y = rows(sol, n, m);

    // Set the statistics of the calculation
    result.statistics.converged = true;
}

auto KktSolver::Impl::solveFullspaceSparse(const KktProblem& problem, KktResult& result, const KktOptions& options) -> void
{
    error("Cannot solve the KKT problem.", "The FullspaceSparse algorithm has not been implemented yet.");

    result.statistics.converged = true;
}

auto KktSolver::Impl::solveRangespace(const KktProblem& problem, KktResult& result, const KktOptions& options) -> void
{
    // The references to the components of the KKT problem
    const auto& A = has_any_scaling ? scaled_problem.A : problem.A;
    const auto& H = has_any_scaling ? scaled_problem.H : problem.H;
    const auto& f = has_any_scaling ? scaled_problem.f : problem.f;
    const auto& g = has_any_scaling ? scaled_problem.g : problem.g;

    // The references to the solution of the KKT problem
    auto& x = result.solution.x;
    auto& y = result.solution.y;

    if(options.diagonalH)
    {
        invH_diag = H.diagonal().cwiseInverse();
        AinvH = A * diag(invH_diag);
        AinvHAt = AinvH * A.transpose();

        llt_AinvHAt.compute(AinvHAt);

        if(llt_AinvHAt.info() != Eigen::Success)
            error("Cannot solve the KKT problem with the rangespace algorithm.",
                "The provided matrix H might not be symmetric positive-definite.");

        y = llt_AinvHAt.solve(g - AinvH*f);
        x = diag(invH_diag) * f + AinvH.transpose()*y;
    }
    else
    {
        error("Cannot solve the KKT problem.", "The Rangespace algorithm has been implemented only for diagonal matrices H.");
    }

    result.statistics.converged = true;
}

auto KktSolver::Impl::solveNullspace(const KktProblem& problem, KktResult& result, const KktOptions& options) -> void
{
    // Update the nullspace and rangespace matrices Z and Y
    computeZY(problem, options);

    // The references to the components of the KKT problem
    const auto& H = has_any_scaling ? scaled_problem.H : problem.H;
    const auto& f = has_any_scaling ? scaled_problem.f : problem.f;
    const auto& g = has_any_scaling ? scaled_problem.g : problem.g;
    const auto& Z = has_any_scaling ? scaled_Z : this->Z;
    const auto& Y = has_any_scaling ? scaled_Y : this->Y;

    // The references to the solution of the KKT problem
    auto& x = result.solution.x;
    auto& y = result.solution.y;

    const Matrix Heff = 0.5 * (H + H.transpose());

    Matrix ZTHZ = Z.transpose() * Heff * Z;

    LLT<Matrix> llt_ZTHZ(ZTHZ);

    if(llt_ZTHZ.info() != Eigen::Success)
        error("Cannot solve the KKT equation using the nullspace algorithm.",
            "The provided H matrix might not be symmetric positive-definite or is ill-conditioned.");

    Vector tmp = Z.transpose() * (f - Heff*Y*g);

    Vector xZ = llt_ZTHZ.solve(tmp);

    x = Z*xZ + Y*g;
    y = Y.transpose() * (Heff*x - f);

    result.statistics.converged = true;
}

auto KktSolver::Impl::solveNullspacePartial(const KktProblem& problem, KktResult& result, const KktOptions& options) -> void
{
    error("Cannot solve the KKT problem.", "The NullspacePartial algorithm has not been implemented yet.");

    result.statistics.converged = true;
}

auto KktSolver::Impl::solve(const KktProblem& problem, KktResult& result, const KktOptions& options) -> void
{
    Time begin = time();

    initialise(problem, options);

    switch(options.algorithm)
    {
    case KktFullspaceDense    : solveFullspaceDense(problem, result, options); break;
    case KktFullspaceSparse   : solveFullspaceSparse(problem, result, options); break;
    case KktRangespace        : solveRangespace(problem, result, options); break;
    case KktNullspace         : solveNullspace(problem, result, options); break;
    case KktNullspacePartial  : solveNullspacePartial(problem, result, options); break;
    default: error("Cannot solve the KKT problem.", "Unknown algorithm.");
    }

    // Compute the unscaled solution vectors `x` and `y`
    if(has_xscaling) result.solution.x = options.xscaling % result.solution.x;
    if(has_yscaling) result.solution.y = options.yscaling % result.solution.y;

    result.statistics.time = elapsed(begin);
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

auto KktSolver::solve(const KktProblem& problem, KktResult& result) -> void
{
    KktOptions default_options;
    pimpl->solve(problem, result, default_options);
}

auto KktSolver::solve(const KktProblem& problem, KktResult& result, const KktOptions& options) -> void
{
    pimpl->solve(problem, result, options);
}
} // namespace Reaktor
