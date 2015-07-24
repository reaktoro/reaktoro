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

struct OptimumSolverIpAction::Impl
{
    KktVector rhs;
    KktSolution sol;
    KktSolver kkt;
    ObjectiveResult f;

    Outputter outputter;

    auto solve(OptimumProblem problem, OptimumState& state, OptimumOptions options) -> OptimumResult;

    auto solveMain(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult;
};

auto OptimumSolverIpAction::Impl::solve(OptimumProblem problem, OptimumState& state, OptimumOptions options) -> OptimumResult
{
    // The transpose of the coefficient matrix `A`
    const Matrix At = tr(problem.A);

    // Calculate the QR decomposition of the transpose of `A`
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(At);

    // Identify the indices of the linearly independent rows of `A`
    const unsigned rank = qr.rank();
    Eigen::VectorXi I = qr.colsPermutation().indices().segment(0, rank);
    std::sort(I.data(), I.data() + rank);

    // The indices of the linearly independent rows of `A`
    const Indices ic(I.data(), I.data() + rank);

    // Define the regularized optimization problem without linearly dependent constraints
    problem.A = rows(problem.A, ic);
    problem.b = rows(problem.b, ic);

    // Remove the names of the linearly dependent constraints
    if(options.output.ynames.size())
        options.output.ynames = extract(options.output.ynames, ic);

    // Get the linearly independent components of the Lagrange multipliers `y`
    state.y = rows(state.y, ic);

    // Solve the regularized optimization problem
    auto result = solveMain(problem, state, options);

    // Calculate the Lagrange multipliers for all equality constraints
    state.y = qr.solve(f.grad - state.z);

    return result;
}

auto OptimumSolverIpAction::Impl::solveMain(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
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
    const auto& tau_x = options.ipaction.tau_x;
    const auto& tau_z = options.ipaction.tau_z;

    // Define some auxiliary references to variables
    auto& x = state.x;
    auto& y = state.y;
    auto& z = state.z;

    // Calculate the LU factorization of the coefficient matrix A
    auto lu = problem.A.fullPivLu();

    // Calculate the kernel (nullspace) matrix K of A such that A*K = 0
    const Matrix K = lu.kernel();

    // Get the lower and upper matrices
    Matrix L = lu.matrixLU().leftCols(m).triangularView<Eigen::UnitLower>();
    Matrix U = lu.matrixLU().triangularView<Eigen::Upper>();

    // Get the permutation matrices
    const auto P1 = lu.permutationP();
    const auto P2 = lu.permutationQ();

    // Set the U1 and U2 submatrices of U = [U1 U2]
    const Matrix U1 = U.leftCols(m);
    const Matrix U2 = U.rightCols(n - m);

    // Compute the regularized coefficient matrix A (leave it as is - do not clean round-off errors)
    const Matrix A = U*P2.inverse();

    // Compute the regularized vector b
    const Vector b = L.triangularView<Eigen::UnitLower>().solve(P1 * problem.b);

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

    // The function that computes the Action step
    auto compute_newton_step = [&]()
    {
        Matrix J = zeros(n, n);
        block(J, 0, 0, n - m, n) = tr(K) * diag(f.hessian.diagonal + z/x);
        block(J, n - m, 0, m, n) = A;

        Matrix r = zeros(n);
        rows(r, 0, n - m) = -tr(K) * (f.grad - mu/x);
        rows(r, n - m, m) = -h;

        sol.dx = J.lu().solve(r);
        sol.dz = (mu - z % x - z % sol.dx)/x;
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
        alphax = fractionToTheBoundary(x, sol.dx, tau_x);
        alphaz = fractionToTheBoundary(z, sol.dz, tau_z);

        x += alphax * sol.dx;
        z += alphaz * sol.dz;
    };

    // The function that computes the current error norms
    auto update_errors = [&]()
    {
        // Calculate the optimality, feasibility and centrality errors
        errorf = norminf(tr(K) * (f.grad - z));
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

    outputter.outputHeader();

    // Finish timing the calculation
    result.time = elapsed(begin);

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
