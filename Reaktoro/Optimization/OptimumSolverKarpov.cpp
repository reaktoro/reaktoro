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

#include "OptimumSolverKarpov.hpp"

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

auto largestStepSize(const Vector& x, const Vector& dx) -> double
{
    double alpha = infinity();
    for(unsigned i = 0; i < x.size(); ++i)
        if(dx[i] < 0.0) alpha = std::min(alpha, -x[i]/dx[i]);
    return alpha;
}

} // namespace

struct OptimumSolverKarpov::Impl
{
    // The result of the objective function evaluation
    ObjectiveResult f;

    // The vector defined as `t = tr(A)*y - grad(f)`
    Vector t;

    // The descent direction vector
    Vector dx;

    // The left-hand and right-hand side matrix and vector of the linear system
    Matrix lhs;
    Vector rhs;

    // The outputter instance
    Outputter outputter;

    auto solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult;
};

auto OptimumSolverKarpov::Impl::solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    // Start timing the calculation
    Time begin = time();

    // Initialize the outputter instance
    outputter = Outputter();
    outputter.setOptions(options.output);

    // The result of the calculation
    OptimumResult result;

    // The number of primal variables `n` and equality constraints `m`
    const unsigned n = problem.A.cols();
    const unsigned m = problem.A.rows();

    // Auxiliary references to primal and dual variables
    Vector& x = state.x;
    Vector& y = state.y;
    Vector& z = state.z;

    // The Lagrange multiplier with respect to the ellipsoid constraint
    double p = 0;

    // Ensure the dual variables `y` and `z` are initialized
    if(y.rows() == 0) y = zeros(m);
    if(z.rows() == 0) z = zeros(n);

    // Auxiliary references to the components of the equality linear constraints
    const auto& A = problem.A;
    const auto& b = problem.b;

    // Auxiliary references to algorithm parameters
    const auto& tolerance = options.tolerance;
    const auto& max_iterations = options.max_iterations;
    const auto& line_search_algorithm = options.karpov.line_search_algorithm;
    const auto& line_search_tolerance = options.karpov.line_search_tolerance;
    const auto& line_search_max_iterations = options.karpov.line_search_max_iterations;
    const auto& line_search_upper_bound = options.karpov.line_search_upper_bound;
    const auto& line_search_upper_factor = options.karpov.line_search_upper_factor;
    const auto& feasibility_tolerance = options.karpov.feasibility_tolerance;
    const auto& tau_feasible = options.karpov.tau_feasible;
    const auto& tau_descent = options.karpov.tau_descent;

    // The alpha step sizes used to restric the steps inside the feasible domain
    double alpha_max, alpha;

    // The current residual error of the calculation
    double error;

    // The current infeasibility error of the calculation
    double infeasibility;

    // The current iteration number
    unsigned iter = 1;

    // The function that outputs the header and initial state of the solution
    auto output_header = [&]()
    {
        if(not options.output.active) return;

        outputter.addEntry("iter");
        outputter.addEntries(options.output.xprefix, n, options.output.xnames);
        outputter.addEntries(options.output.yprefix, m, options.output.ynames);
        outputter.addEntries(options.output.zprefix, n, options.output.znames);
        outputter.addEntry("p");
        outputter.addEntry("f(x)");
        outputter.addEntry("h(x)");
        outputter.addEntry("error");
        outputter.addEntry("alpha");
        outputter.addEntry("alpha[upper]");
        outputter.outputHeader();

        // Evaluate the objective function at the initial guess `x`
        f = problem.objective(x);

        // Calculate the initial infeasibility
        infeasibility = norm(A*x - b);

        outputter.addValue(result.iterations);
        outputter.addValues(x);
        outputter.addValues(y);
        outputter.addValues(z);
        outputter.addValue(p);
        outputter.addValue(f.val);
        outputter.addValue(infeasibility);
        outputter.addValue("---");
        outputter.addValue("---");
        outputter.addValue("---");
        outputter.outputState();
    };

    // The function that outputs the current state of the solution
    auto output_state = [&]()
    {
        if(not options.output.active) return;

        outputter.addValue(iter);
        outputter.addValues(x);
        outputter.addValues(y);
        outputter.addValues(z);
        outputter.addValue(p);
        outputter.addValue(f.val);
        outputter.addValue(infeasibility);
        outputter.addValue(error);
        outputter.addValue(alpha);
        outputter.addValue(alpha_max);
        outputter.outputState();
    };

    // Initialize the variables before the calculation begins
    auto initialize = [&]()
    {
        // Ensure the line search step length starts at its upper bound
        alpha = line_search_upper_bound;
    };

    // Calculate a feasible point for the minimization calculation
    auto calculate_feasible_point = [&]()
    {
        for(; iter < max_iterations; ++iter)
        {
            // Check if the current iterate is sufficiently feasible
            if(infeasibility < feasibility_tolerance)
                break;

            // Assemble the linear system whose rhs is the feasibility residual
            lhs = A*diag(x)*tr(A);
            rhs = b - A*x;

            // Calculate the dual variables `y`
            y = lhs.fullPivLu().solve(rhs);

            // Calculate the correction step towards a feasible point
            dx = diag(x)*tr(A)*y;

            // Calculate the largest step size
            alpha_max = tau_feasible * largestStepSize(x, dx);

            // Ensure the step size is bounded above by 1
            alpha = std::min(1.0, alpha_max);

            // Update the primal variables
            x += alpha * dx;

            // Calculate the feasibility residual
            infeasibility = norm(A*x - b);

            output_state();
        }

        // Evaluate the objective function at the feasible point `x`
        f = problem.objective(x);
    };

    // Calculate the descent direction
    auto calculate_descent_direction = [&]()
    {
        // Assemble the linear system
        lhs = A*diag(x)*tr(A);
        rhs = A*diag(x)*f.grad;

        // Solve for the dual variables `y`
        y = lhs.lu().solve(rhs);

        // Calculate the auxiliary vector `t = tr(A)*y - grad(f)`
        t = tr(A)*y - f.grad;

        // Calculate the dual variable `p`
        p = tr(t)*diag(x)*t;
        p = std::sqrt(p);

        // Calculate the descent step for the primal variables `x`
        dx = diag(x)*t;
    };

    // Solve the line search minimization problem
    auto solve_line_search_minimization_problem = [&]()
    {
        // Calculate the largest step size that will not leave the feasible domain
        alpha_max = tau_descent * largestStepSize(x, dx);

        // Enforce the maximum step length is not too big with respect to the last step length
        alpha_max = std::min(alpha_max, line_search_upper_factor*alpha);

        // Enforce the maximum step length is not over its upper bound
        alpha_max = std::min(alpha_max, line_search_upper_bound);

        // Define the single variable function for the line search minimization problem
        auto g = [&](double alpha) -> double
        {
            return problem.objective(x + alpha*dx).val;
        };

        // Solve the line search minimization problem
        if(line_search_algorithm == "GoldenSectionSearch")
            alpha = minimizeGoldenSectionSearch(g, 0.0, alpha_max, line_search_tolerance);
        else
            alpha = minimizeBrent(g, 0.0, alpha_max, line_search_tolerance, line_search_max_iterations);
    };

    // Update the current state of the calculation
    auto update_state = [&]()
    {
        // Calculate the new primal iterate using the calculated `alpha` step length
        x += alpha * dx;

        // Calculate the new dual iterate `z`
        z = -t;

        // Evaluate the objective function at the new iterate `x`
        f = problem.objective(x);
    };

    // Update the error norms
    auto update_errors = [&]()
    {
        // Calculate the current error of the minimization calculation
        error = norm(diag(x)*t);
    };

    // Return true if the calculation has converged
    auto converged = [&]()
    {
        if(error < tolerance)
        {
            result.iterations = iter;
            result.succeeded = true;
            return true;
        }
        return false;
    };

    initialize();

    output_header();

    calculate_feasible_point();

    for(; iter < max_iterations; ++iter)
    {
        calculate_descent_direction();

        solve_line_search_minimization_problem();

        update_state();

        update_errors();

        output_state();

        if(converged())
        {
            result.iterations = iter;
            result.succeeded = true;
            break;
        }
    }

    // Output the header of the calculation
    outputter.outputHeader();

    // Finish timing the calculation
    result.time = elapsed(begin);

    return result;
}

OptimumSolverKarpov::OptimumSolverKarpov()
: pimpl(new Impl())
{}

OptimumSolverKarpov::OptimumSolverKarpov(const OptimumSolverKarpov& other)
: pimpl(new Impl(*other.pimpl))
{}

OptimumSolverKarpov::~OptimumSolverKarpov()
{}

auto OptimumSolverKarpov::operator=(OptimumSolverKarpov other) -> OptimumSolverKarpov&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto OptimumSolverKarpov::solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    return pimpl->solve(problem, state, options);
}

} // namespace Reaktoro
