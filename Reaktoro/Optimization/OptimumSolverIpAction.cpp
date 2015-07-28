// Reaktoro is a C++ library for computational reaction modelling.
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

#include "OptimumSolverIpAction.hpp"

// Eigen includes
#include <eigen/Dense>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Outputter.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/TimeUtils.hpp>
#include <Reaktoro/Math/MathUtils.hpp>
#include <Reaktoro/Optimization/KktSolver.hpp>
#include <Reaktoro/Optimization/OptimumProblem.hpp>
#include <Reaktoro/Optimization/OptimumOptions.hpp>
#include <Reaktoro/Optimization/OptimumResult.hpp>
#include <Reaktoro/Optimization/OptimumState.hpp>
#include <Reaktoro/Optimization/Utils.hpp>

namespace Reaktoro {
namespace {

struct LinearSystemDiagonalSolver
{
    Vector A1, A2, invA1, a1, a2, q, u;
    Matrix B1, B2, C1, C2, Q;
    Indices ipivot, inonpivot;

    auto solve(
        const Vector& A,
        const Matrix& B,
        const Matrix& C,
        const Matrix& D,
        const Vector& a,
        const Vector& b,
        Vector& x,
        Vector& y) -> void
    {
        const unsigned n = A.rows();
        const unsigned m = B.cols();
        ipivot.clear();
        inonpivot.clear();
        ipivot.reserve(n);
        inonpivot.reserve(n);
        for(unsigned i = 0; i < n; ++i)
            if(A[i] > norminf(C.col(i))) ipivot.push_back(i);
            else inonpivot.push_back(i);

        const unsigned n2 = inonpivot.size();

        A1 = rows(A, ipivot);
        A2 = rows(A, inonpivot);
        invA1 = 1/A1;
        B1 = rows(B, ipivot);
        B2 = rows(B, inonpivot);
        C1 = cols(C, ipivot);
        C2 = cols(C, inonpivot);
        a1 = rows(a, ipivot);
        a2 = rows(a, inonpivot);

        Q = zeros(n2 + m, n2 + m);
        Q.topLeftCorner(n2, n2).diagonal() = A2;
        Q.topRightCorner(n2, m) = B2;
        Q.bottomLeftCorner(m, n2) = C2;
        Q.bottomRightCorner(m, m).noalias() = D - C1*diag(invA1)*B1;

        q.resize(n2 + m);
        rows(q, 0, n2) = a2;
        rows(q, n2, m) = b - C1*diag(invA1)*a1;

        u = Q.lu().solve(q);

        y = rows(u, n2, m);

        x.resize(n);
        rows(x, ipivot) = invA1 % (a1 - B1*y);
        rows(x, inonpivot) = rows(u, 0, n2);
    }
};

} // namespace

struct OptimumSolverIpAction::Impl
{
    KktVector rhs;
    KktSolution sol;
    KktSolver kkt;

    Outputter outputter;

    /// The dedicated linear system solver with a diagonal matrix on the top-left corner
    LinearSystemDiagonalSolver lsd;

    /// The coefficient matrix `A` from the last calculation
    Matrix lastA;

    /// The LU factorization of the coefficient matrix `A`
    Eigen::FullPivLU<Matrix> lu;

    /// The first `m` indices of the permutation matrix `Q`
    Indices iQ1;

    /// The last `n - m` indices of the permutation matrix `Q`
    Indices iQ2;

    /// The kernel (nullspace) matrix `K` of `tr(A)` such that `K*tr(A) = 0`
    Matrix K;

    /// The matrix formed from the columns of `K` that corresponds to the indices in `iQ1`
    Matrix K1;

    /// The regularized coefficient matrix `A`
    Matrix A;

    /// The regularized right-hand side vector `b`
    Vector b;

    // Update the matrices `A1` and `A2` whose columns correspond to the indices `iQ1` and `iQ2`
    Matrix A1, A2;

    auto solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult;
};

auto OptimumSolverIpAction::Impl::solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    // Start timing the calculation
    Time begin = time();

    // Initialize the outputter instance
    outputter = Outputter();
    outputter.setOptions(options.output);

    // Set the KKT options
    kkt.setOptions(options.kkt);

    // The result of the calculation
    OptimumResult result;

    // The number of primal variables and equality constraints
    const unsigned n = problem.A.cols();
    const unsigned m = problem.A.rows();

    // Define some auxiliary references to parameters
    const auto& tolerance = options.tolerance;
    const auto& mu = options.ipaction.mu;
    const auto& tau = options.ipaction.tau;

    // Define some auxiliary references to variables
    auto& x = state.x;
    auto& y = state.y;
    auto& z = state.z;
    auto& f = state.f;

    // Ensure the initial guesses for `x` and `y` have adequate dimensions
    if(x.size() != n) x = zeros(n);
    if(y.size() != m) y = zeros(m);
    if(z.size() != n) z = zeros(n);

    // Ensure the initial guesses for `x` and `z` are inside the feasible domain
    x = (x.array() > 0.0).select(x, 1.0);
    z = (z.array() > 0.0).select(z, 1.0);

    // The auxiliary vector containing the feasibility residual Ax - b
    Vector h;

    // The alpha step sizes used to restric the steps inside the feasible domain
    double alphax, alphaz;

    // The optimality, feasibility, centrality and total error variables
    double errorf, errorh, errorc, error;

    // The function that outputs the header and initial state of the solution
    auto output_header = [&]()
    {
        if(!options.output.active) return;

        outputter.addEntry("iter");
        outputter.addEntries(options.output.xprefix, n, options.output.xnames);
        outputter.addEntries(options.output.zprefix, n, options.output.znames);
        outputter.addEntry("f(x)");
        outputter.addEntry("h(x)");
        outputter.addEntry("errorf");
        outputter.addEntry("errorh");
        outputter.addEntry("errorc");
        outputter.addEntry("error");
        outputter.addEntry("alphax");
        outputter.addEntry("alphaz");

        outputter.outputHeader();
        outputter.addValue(result.iterations);
        outputter.addValues(x);
        outputter.addValues(z);
        outputter.addValue(f.val);
        outputter.addValue(norminf(h));
        outputter.addValue("---");
        outputter.addValue("---");
        outputter.addValue("---");
        outputter.addValue("---");
        outputter.addValue("---");
        outputter.addValue("---");
        outputter.outputState();
    };

    // The function that outputs the current state of the solution
    auto output_state = [&]()
    {
        if(!options.output.active) return;

        outputter.addValue(result.iterations);
        outputter.addValues(x);
        outputter.addValues(z);
        outputter.addValue(f.val);
        outputter.addValue(norminf(h));
        outputter.addValue(errorf);
        outputter.addValue(errorh);
        outputter.addValue(errorc);
        outputter.addValue(error);
        outputter.addValue(alphax);
        outputter.addValue(alphaz);
        outputter.outputState();
    };

    // The function that initialize the solver
    auto initialize = [&]()
    {
        // Skip the update steps below if matrix `A` has not changed
        if(lastA.rows() != m || lastA.cols() != n || problem.A != lastA)
        {
            // Update the last given coefficient matrix `A`
            lastA = A;

            // Update the LU factorization of the coefficient matrix A
            lu.compute(problem.A);

            // Get the lower and upper matrices
            const Matrix U = lu.matrixLU().triangularView<Eigen::Upper>();

            // Get the permutation matrix Q, where PAQ = LU
            const auto Q = lu.permutationQ();

            // Update the first `m` indices of the permutation matrix `Q`
            iQ1 = Indices(Q.indices().data(), Q.indices().data() + m);

            // Update the last `n - m` indices of the permutation matrix `Q`
            iQ2 = Indices(Q.indices().data() + m, Q.indices().data() + n);

            // Update the kernel (nullspace) matrix K of tr(A) such that K*tr(A) = 0
            K = tr(lu.kernel());

            // Update the regularized coefficient matrix A (leave it as is - do not clean round-off errors)
            A = U * Q.inverse();

            // Update the matrices `A1` and `A2` whose columns correspond to the indices `iQ1` and `iQ2`
            A1 = cols(A, iQ1);
            A2 = cols(A, iQ2);

            // Update the matrix `K1` formed from the colums of the kernel matrix `K` corresponding to the indices `iQ1`
            K1 = cols(K, iQ1);
        }

        // Get the permutation matrix `P`, where `PAQ = LU`
        const auto P = lu.permutationP();

        // Get the lower factor of the coefficient matrix `A`
        const auto L = lu.matrixLU().leftCols(m).triangularView<Eigen::UnitLower>();

        // Update the regularized vector b
        b = L.solve(P * problem.b);
    };

    // The function that updates the objective and constraint state
    auto update_state = [&]()
    {
        f = problem.objective(x);
        h = A*x - b;
    };

    // Return true if function `update_state` failed
    auto update_state_failed = [&]()
    {
        const bool f_finite = std::isfinite(f.val);
        const bool g_finite = f.grad.allFinite();
        const bool all_finite = f_finite && g_finite;
        return !all_finite;
    };

    // The function that computes the Newton step
    auto compute_newton_step_diagonal = [&]()
    {
        Vector D = f.hessian.diagonal + z/x;
        Vector D1 = rows(D, iQ1);
        Vector D2 = rows(D, iQ2);

        Matrix B = K1 * diag(D1);

        Vector dx1, dx2;

        const Vector m1 = -K*(f.grad - mu/x);
        const Vector m2 = -(A*x - b);

        lsd.solve(D2, B, A2, A1, m1, m2, dx2, dx1);

        sol.dx.resize(n);
        rows(sol.dx, iQ1) = dx1;
        rows(sol.dx, iQ2) = dx2;

        sol.dz = (mu - z % x - z % sol.dx)/x;
    };

    // The function that computes the Newton step
    auto compute_newton_step_dense = [&]()
    {
        f.hessian.dense.diagonal() += z/x;

        Matrix J = zeros(n, n);
        block(J, 0, 0, n - m, n) = K*f.hessian.dense;
        block(J, n - m, 0, m, n) = A;

        Vector r = zeros(n);
        rows(r, 0, n - m) = -K*(f.grad - mu/x);
        rows(r, n - m, m) = -h;

        sol.dx = J.lu().solve(r);
        sol.dz = (mu - z % x - z % sol.dx)/x;
    };

    auto compute_newton_step = [&]()
    {
        if(f.hessian.mode == Hessian::Diagonal)
            compute_newton_step_diagonal();
        else
            compute_newton_step_dense();
    };

    // Return true if the function `compute_newton_step` failed
    auto compute_newton_step_failed = [&]()
    {
        const bool dx_finite = sol.dx.allFinite();
        const bool dz_finite = sol.dz.allFinite();
        const bool all_finite = dx_finite && dz_finite;
        return !all_finite;
    };

    // The function that performs an update in the iterates
    auto update_iterates = [&]()
    {
        alphax = fractionToTheBoundary(x, sol.dx, tau);
        alphaz = fractionToTheBoundary(z, sol.dz, tau);

        x += alphax * sol.dx;
        z += alphaz * sol.dz;
    };

    // The function that computes the current error norms
    auto update_errors = [&]()
    {
        // Calculate the optimality, feasibility and centrality errors
        errorf = norminf(K*(f.grad - z));
        errorh = norminf(h);
        errorc = norminf(x%z - mu);

        // Calculate the maximum error
        error = std::max({errorf, errorh, errorc});
        result.error = error;
    };

    auto converged = [&]()
    {
        if(error < tolerance)
        {
            result.succeeded = true;
            return true;
        }
        return false;
    };

    initialize();
    update_state();
    output_header();

    do
    {
        ++result.iterations; if(result.iterations > options.max_iterations) break;
        compute_newton_step();
        if(compute_newton_step_failed())
            break;
        update_iterates();
        update_state();
        if(update_state_failed())
            break;
        update_errors();
        output_state();
    } while(!converged());

    // Finish timing the calculation
    result.time = elapsed(begin);

    outputter.outputHeader();
    outputter.outputMessage("Calculation completed in ", result.time, " seconds.");

    return result;
}

OptimumSolverIpAction::OptimumSolverIpAction()
: pimpl(new Impl())
{}

OptimumSolverIpAction::OptimumSolverIpAction(const OptimumSolverIpAction& other)
: pimpl(new Impl(*other.pimpl))
{}

OptimumSolverIpAction::~OptimumSolverIpAction()
{}

auto OptimumSolverIpAction::operator=(OptimumSolverIpAction other) -> OptimumSolverIpAction&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto OptimumSolverIpAction::solve(const OptimumProblem& problem, OptimumState& state) -> OptimumResult
{
    return pimpl->solve(problem, state, {});
}

auto OptimumSolverIpAction::solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    return pimpl->solve(problem, state, options);
}

} // namespace Reaktoro
