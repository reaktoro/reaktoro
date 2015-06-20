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

#include "OptimumSolverIpActive.hpp"

// Eigen includes
#include <eigen/Dense>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Outputter.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/TimeUtils.hpp>
#include <Reaktoro/Math/MathUtils.hpp>
#include <Reaktoro/Optimization/KktSolver.hpp>
#include <Reaktoro/Optimization/OptimumOptions.hpp>
#include <Reaktoro/Optimization/OptimumProblem.hpp>
#include <Reaktoro/Optimization/OptimumResult.hpp>
#include <Reaktoro/Optimization/OptimumSolverIpNewton.hpp>
#include <Reaktoro/Optimization/OptimumState.hpp>
#include <Reaktoro/Optimization/Utils.hpp>

namespace Reaktoro {
namespace {

auto throwZeroInitialGuessError() -> void
{
    Exception exception;
    exception.error << "Cannot continue the optimization calculation with given initial guess.";
    exception.reason << "The provided initial guess has all primal variables on their bounds.";
    RaiseError(exception);
}

} // namespace

struct OptimumSolverIpActive::Impl
{
    KktVector rhs;
    KktSolution sol;
    KktSolver kkt;
    ObjectiveResult f;

    /// The solver for the optimisation calculations
    OptimumSolverIpNewton solver;

    /// The stable primal variables
    Vector xs;

    /// The dual potentials of the unstable primal variables
    Vector zu;

    /// The gradient of the objective function w.r.t. stable primal variables
    Vector gs;

    /// The gradient of the objective function w.r.t. stable primal variables
    Vector gu;

    /// The indices of the stable primal variables
    Indices istable_variables;

    /// The indices of the unstable primal variables
    Indices iunstable_variables;

    /// The coefficient matrix of the stable variables
    Matrix As;

    /// The coefficient matrix of the unstable variables
    Matrix Au;

    /// The optimisation problem based on stable variables only
    OptimumProblem stable_problem;

    /// The state of the optimisation calculation based on stable variables only
    OptimumState stable_state;

    /// The options for the optimisation calculation based on stable variables only
    OptimumOptions stable_options;

    // The outputter instance
    Outputter outputter;

    auto solve(OptimumProblem problem, OptimumState& state, OptimumOptions options) -> OptimumResult;

    auto solveMain(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult;
};

auto OptimumSolverIpActive::Impl::solve(OptimumProblem problem, OptimumState& state, OptimumOptions options) -> OptimumResult
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

auto OptimumSolverIpActive::Impl::solveMain(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
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
    const auto& n = problem.A.cols();
    const auto& m = problem.A.rows();
    const auto& zero = options.ipactive.epsilon;

    // Ensure the initial guesses for `x` and `y` have adequate dimensions
    if(x.size() != n) x = zeros(n);
    if(y.size() != m) y = zeros(m);
    if(z.size() != n) z = zeros(n);

    // Ensure the initial guesses for `x` and `z` are inside the feasible domain
    x = (x.array() > 0.0).select(x, 1.0);
    z = (z.array() > 0.0).select(z, 1.0);

    // Update the sets of stable and unstable variables
    // by checking which variables have zero molar amounts
    auto update_stability_sets = [&]()
    {
        // Update the indices of the stable and unstable variables
        istable_variables.clear();
        iunstable_variables.clear();
        for(unsigned i = 0; i < n; ++i)
            if(x[i] > 0.0) istable_variables.push_back(i);
            else iunstable_variables.push_back(i);

        // Check if the set of stable variables is empty
        if(istable_variables.empty())
            throwZeroInitialGuessError();

        // Update the balance matrices of stable and unstable variables
        As = cols(A, istable_variables);
        Au = cols(A, iunstable_variables);
    };

    /// Update the stable primal variables by setting them to zero
    auto update_unstable_variables = [&]()
    {
        for(auto i : istable_variables)
            if(x[i] < zero)
                x[i] = 0.0;
    };

    /// Update the OptimumOptions instance `optimum_options`
    auto update_optimum_options = [&]()
    {
        // Initialize the options for the optimisation calculation
        stable_options = options;
        stable_options.ipnewton.mu = zero * 1e-5;

        // Initialize the names of the primal and dual variables
        if(options.output.active)
        {
            // Create references to the primal and dual variables `x` and `z`
            auto& xnames = stable_options.output.xnames;
            auto& znames = stable_options.output.znames;

            // Ensure only the names of stable primal variables are output
            if(xnames.size())
            {
                xnames = extract(xnames, istable_variables);
                znames = xnames;
            }
        }
    };

    /// Update the OptimumProblem instance `optimum_problem`
    auto update_optimum_problem = [&]()
    {
        // The number of stable variables and elements in the equilibrium partition
        const unsigned num_stable_variables = istable_variables.size();

        // The result of the objective evaluation
        ObjectiveResult res;

        stable_problem.objective = [=](const Vector& xs) mutable
        {
            // Update the stable components in `x`
            rows(x, istable_variables) = xs;

            // Evaluate the objective function using updated `x`
            res = problem.objective(x);

            return res;
        };

        stable_problem.A = As;
        stable_problem.b = b;
        stable_problem.l = zeros(num_stable_variables);
    };

    /// Update the OptimumState instance `stable_state` using `state`
    auto update_stable_state = [&]()
    {
        // Update the optimum state of the stable components
        rows(x, istable_variables).to(stable_state.x);
        rows(z, istable_variables).to(stable_state.z);
        stable_state.y = y;
    };

    /// Update the OptimumState instance `state` using `stable_state`
    auto update_state = [&]()
    {
        // Set all dual variables `z` to zero, but correct the stable components below
        z.fill(0.0);

        // Update the stable components of the primal and dual variables `x` and `z`
        rows(x, istable_variables) = stable_state.x;
        rows(z, istable_variables) = stable_state.z;

        // Update the dual variables `y` with respect to the linear constraints
        y = stable_state.y;
    };

    /// Return true if not all variables are still stable
    auto not_all_stable_yet = [&]()
    {
        if(iunstable_variables.empty())
            return true;

        const double lambda = std::sqrt(zero);

        for(Index i : iunstable_variables)
            x[i] = zero;

        f = problem.objective(x);

        gu = rows(f.grad, iunstable_variables);

        zu = gu - tr(Au) * y;

        for(int k = 0; k < zu.rows(); ++k)
            if(zu[k] < 0.0)
                state.x[iunstable_variables[k]] = lambda;

        return min(zu) < 0.0;
    };

    do {
        // Update the sets of stable and unstable variables
        update_stability_sets();

        // Update the optimum options
        update_optimum_options();

        // Update the optimum problem
        update_optimum_problem();

        // Update the optimum state
        update_stable_state();

        // Solve the optimisation problem
        result += solver.solve(stable_problem, stable_state, stable_options);

        // Exit if the the optimisation calculation did not succeed
        if(!result.succeeded) break;

        // Update the full state from the stable state
        update_state();

    } while(not_all_stable_yet());

    update_unstable_variables();

    // Finish timing the calculation
    result.time = elapsed(begin);

    return result;
}

OptimumSolverIpActive::OptimumSolverIpActive()
: pimpl(new Impl())
{}

OptimumSolverIpActive::OptimumSolverIpActive(const OptimumSolverIpActive& other)
: pimpl(new Impl(*other.pimpl))
{}

OptimumSolverIpActive::~OptimumSolverIpActive()
{}

auto OptimumSolverIpActive::operator=(OptimumSolverIpActive other) -> OptimumSolverIpActive&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto OptimumSolverIpActive::solve(const OptimumProblem& problem, OptimumState& state) -> OptimumResult
{
    return pimpl->solve(problem, state, {});
}

auto OptimumSolverIpActive::solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    return pimpl->solve(problem, state, options);
}

} // namespace Reaktoro
