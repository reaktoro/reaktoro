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

#include "OptimumSolverIpnewton.hpp"

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

struct OptimumSolverIpnewton::Impl
{
    KktVector rhs;
    KktSolution sol;
    KktSolver kkt;
    ObjectiveResult f;

    Outputter outputter;

    auto solve(OptimumProblem problem, OptimumState& state, OptimumOptions options) -> OptimumResult;

    auto solveMain(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult;
};

auto OptimumSolverIpnewton::Impl::solve(OptimumProblem problem, OptimumState& state, OptimumOptions options) -> OptimumResult
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

auto OptimumSolverIpnewton::Impl::solveMain(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
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

    // Define some auxiliary references to variables
    auto& x = state.x;
    auto& y = state.y;
    auto& z = state.z;

    const auto& A = problem.A;
    const auto& b = problem.b;

    Vector h;

    // Define some auxiliary references to parameters
    const auto& n         = problem.A.cols();
    const auto& m         = problem.A.rows();
    const auto& tolerance = options.tolerance;
    const auto& mu        = options.ipnewton.mu;
    const auto& tau       = options.ipnewton.tau;

    // Ensure the initial guesses for `x` and `y` have adequate dimensions
    if(x.size() != n) x = zeros(n);
    if(y.size() != m) y = zeros(m);
    if(z.size() != n) z = zeros(n);

    // Ensure the initial guesses for `x` and `z` are inside the feasible domain
    x = (x.array() > 0.0).select(x, 1.0);
    z = (z.array() > 0.0).select(z, 1.0);

    // The transpose representation of matrix `A`
    const auto At = tr(A);

    // The alpha step sizes used to restric the steps inside the feasible domain
    double alphax, alphaz, alpha;

    // The optimality, feasibility, centrality and total error variables
    double errorf, errorh, errorc, error;

    // The function that outputs the header and initial state of the solution
    auto output_header = [&]()
    {
        if(not options.output.active) return;

        outputter.addEntry("iter");
        outputter.addEntries(options.output.xprefix, n, options.output.xnames);
        outputter.addEntries(options.output.yprefix, m, options.output.ynames);
        outputter.addEntries(options.output.zprefix, n, options.output.znames);
        outputter.addEntry("f(x)");
        outputter.addEntry("h(x)");
        outputter.addEntry("errorf");
        outputter.addEntry("errorh");
        outputter.addEntry("errorc");
        outputter.addEntry("error");
        outputter.addEntry("alpha");
        outputter.addEntry("alphax");
        outputter.addEntry("alphaz");

        outputter.outputHeader();
        outputter.addValue(result.iterations);
        outputter.addValues(x);
        outputter.addValues(y);
        outputter.addValues(z);
        outputter.addValue(f.val);
        outputter.addValue(norminf(h));
        outputter.addValue("---");
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
        if(not options.output.active) return;

        outputter.addValue(result.iterations);
        outputter.addValues(x);
        outputter.addValues(y);
        outputter.addValues(z);
        outputter.addValue(f.val);
        outputter.addValue(norminf(h));
        outputter.addValue(errorf);
        outputter.addValue(errorh);
        outputter.addValue(errorc);
        outputter.addValue(error);
        outputter.addValue(alpha);
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
        const bool all_finite = f_finite and g_finite;
        return not all_finite;
    };

    // The function that computes the Newton step
    auto compute_newton_step = [&]()
    {
        KktMatrix lhs{f.hessian, A, x, z};

        kkt.decompose(lhs);

        // Compute the right-hand side vectors of the KKT equation
        rhs.rx.noalias() = -(f.grad - At*y - z);
        rhs.ry.noalias() = -(h);
        rhs.rz.noalias() = -(x % z - mu);

        // Compute `dx` and `dy` by solving the KKT equation
        kkt.solve(rhs, sol);

        // Update the time spent in linear systems
        result.time_linear_systems += kkt.result().time_solve;
        result.time_linear_systems += kkt.result().time_decompose;
    };

    // Return true if the function `compute_newton_step` failed
    auto compute_newton_step_failed = [&]()
    {
        const bool dx_finite = sol.dx.allFinite();
        const bool dy_finite = sol.dy.allFinite();
        const bool dz_finite = sol.dz.allFinite();
        const bool all_finite = dx_finite and dy_finite and dz_finite;
        return not all_finite;
    };

    // The function that performs an update in the iterates
    auto update_iterates = [&]()
    {
        alphax = fractionToTheBoundary(x, sol.dx, tau);
        alphaz = fractionToTheBoundary(z, sol.dz, tau);
        alpha  = std::min(alphax, alphaz);

        if(options.ipnewton.uniform_newton_step)
        {
            x += alpha * sol.dx;
            y += alpha * sol.dy;
            z += alpha * sol.dz;
        }
        else
        {
            x += alphax * sol.dx;
            y += sol.dy;
            z += alphaz * sol.dz;
        }
    };

    // The function that computes the current error norms
    auto update_errors = [&]()
    {
        // Calculate the optimality, feasibility and centrality errors
        errorf = norminf(f.grad - At*y - z);
        errorh = norminf(h);
        errorc = norminf(x%z/mu - 1);

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
    } while(not converged());

    outputter.outputHeader();

    // Finish timing the calculation
    result.time = elapsed(begin);

    return result;
}

OptimumSolverIpnewton::OptimumSolverIpnewton()
: pimpl(new Impl())
{}

OptimumSolverIpnewton::OptimumSolverIpnewton(const OptimumSolverIpnewton& other)
: pimpl(new Impl(*other.pimpl))
{}

OptimumSolverIpnewton::~OptimumSolverIpnewton()
{}

auto OptimumSolverIpnewton::operator=(OptimumSolverIpnewton other) -> OptimumSolverIpnewton&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto OptimumSolverIpnewton::solve(const OptimumProblem& problem, OptimumState& state) -> OptimumResult
{
    return pimpl->solve(problem, state, {});
}

auto OptimumSolverIpnewton::solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    return pimpl->solve(problem, state, options);
}

} // namespace Reaktoro
