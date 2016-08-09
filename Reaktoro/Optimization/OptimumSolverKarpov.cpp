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

#include "OptimumSolverKarpov.hpp"

// Eigen includes
#include <Reaktoro/Eigen/Cholesky>

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

auto largestStepSize(const Vector& x, const Vector& dx, const Vector& l) -> double
{
    double alpha = infinity();
    for(unsigned i = 0; i < x.size(); ++i)
        if(dx[i] < 0.0)
            alpha = std::min(alpha, std::abs(-(x[i] - l[i])/dx[i]));
    return alpha;
}

} // namespace

struct OptimumSolverKarpov::Impl
{
    // The result of the objective function evaluation
    ObjectiveResult f;

    // The result of the objective function evaluation at a trial step length
    ObjectiveResult f_alpha;

    // The result of the objective function evaluation at the maximum allowed step length
    ObjectiveResult f_alpha_max;

    // The iterate at the trial step length
    Vector x_alpha;

    // The vector of weights for the primal variables in the ellipsoid condition
    Vector w;

    // The vector defined as `t = tr(A)*y - grad(f)`
    Vector t;

    // The descent direction vector
    Vector dx;

    // The left-hand and right-hand side matrix and vector of the linear system
    Matrix lhs;
    Vector rhs;

    // The outputter instance
    Outputter outputter;

    /// Solve the optimization problem
    auto solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
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
        auto& x = state.x;
        auto& y = state.y;
        auto& z = state.z;

        // Auxiliary references to the components of the equality linear constraints
        const auto& A = problem.A;
        const auto& b = problem.b;
        const auto& l = problem.l;

        // The Lagrange multiplier with respect to the ellipsoid constraint
        double p = 0.0;

        // Ensure the primal variables `x` are within the interior domain
        x = max(x, l);

        // Ensure the dual variables `y` and `z` are initialized
        if(y.rows() == 0) y = zeros(m);
        if(z.rows() == 0) z = zeros(n);

        // Auxiliary references to algorithm parameters
        const auto& tolerance = options.tolerance;
        const auto& max_iterations = options.max_iterations;
        const auto& line_search_max_iterations = options.karpov.line_search_max_iterations;
        const auto& line_search_wolfe = options.karpov.line_search_wolfe;
        const auto& feasibility_tolerance = options.karpov.feasibility_tolerance;
        const auto& negative_dual_tolerance = options.karpov.negative_dual_tolerance;
        const auto& active_to_inactive = options.karpov.active_to_inactive;
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
            if(!options.output.active) return;

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
            if(!options.output.active) return;

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
            // Evaluate the objective function at the initial guess `x`
            f = problem.objective(x);

            // Calculate the initial infeasibility
            infeasibility = norm(A*x - b);
        };

        // Calculate a feasible point for the minimization calculation
        auto calculate_feasible_point = [&]()
        {
            outputter.outputMessage("...solving the feasible problem", '\n');

            for(; iter < max_iterations; ++iter)
            {
                // Check if the current iterate is sufficiently feasible
                if(infeasibility < feasibility_tolerance)
                    break;

                // Calculate the weights for the primal variables in the ellipsoid condition
                w.noalias() = x - l;

                // Assemble the linear system whose rhs is the feasibility residual
                lhs = A*diag(w)*tr(A);
                rhs = b - A*x;

                // Calculate the dual variables `y`
                y = lhs.llt().solve(rhs);

                // Calculate the correction step towards a feasible point
                dx = diag(w)*tr(A)*y;

                // Calculate the largest step size
                alpha_max = tau_feasible * largestStepSize(x, dx, l);

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

            outputter.outputMessage("...finished the feasible problem", '\n');
        };

        // Calculate the descent direction using KktSolver
        auto calculate_descent_direction_using_kkt = [&]()
        {
            // Calculate the weights for the primal variables in the ellipsoid condition
            w.noalias() = x - l;

            // Set the Hessian to diagonal mode and use the weights to solve the KKT equation
            f.hessian.mode = Hessian::Diagonal;
            f.hessian.diagonal = inv(w);

            // Set the dual variables `z` to zero
            z = zeros(n);

            // Define the coefficient matrix of the KKT equation
            KktMatrix lhs{f.hessian, A, x, z};

            KktVector rhs;
            rhs.rx = -f.grad;
            rhs.ry = zeros(m);
            rhs.rz = zeros(n);

            KktSolver kkt;
            kkt.setOptions(options.kkt);

            KktSolution sol;

            kkt.decompose(lhs);
            kkt.solve(rhs, sol);

            dx = sol.dx;
            y = sol.dy;
            z = sol.dz;

            // Calculate the auxiliary vector `t = tr(A)*y - grad(f)`
            t = tr(A)*y - f.grad;
        };

        auto calculate_descent_direction_using_llt = [&]()
        {
            // Calculate the weights for the primal variables in the ellipsoid condition
            w.noalias() = x - l;

            // Assemble the linear system
            lhs = A*diag(w)*tr(A);
            rhs = A*diag(w)*f.grad;

            // Solve for the dual variables `y`
            y = lhs.llt().solve(rhs);

            // Calculate the auxiliary vector `t = tr(A)*y - grad(f)`
            t = tr(A)*y - f.grad;

            // Calculate the dual variable `p`
            p = tr(t)*diag(w)*t;
            p = std::sqrt(p);

            // Calculate the descent step for the primal variables `x`
            dx = diag(w)*t;
        };

        auto calculate_descent_direction = [&]()
        {
            if(options.karpov.use_kkt_solver)
                calculate_descent_direction_using_kkt();
            else
                calculate_descent_direction_using_llt();
        };

        // Solve the line search minimization problem
        auto solve_line_search_minimization_problem = [&]()
        {
            // Calculate the largest step size that will not leave the feasible domain
            alpha_max = largestStepSize(x, dx, l);

            // Ensure the largest step size is bounded above by 1
            alpha_max = tau_descent * alpha_max;

            // Start the backtracking line search algorithm
            unsigned i = 0;
            alpha = std::min(alpha_max, 1.0);
            x_alpha = x + alpha*dx;
            f_alpha = f_alpha_max = problem.objective(x_alpha);
            for(; i < line_search_max_iterations; ++i)
            {
                if(!std::isfinite(f_alpha.val) || min(x_alpha - l) < 0.0)
                {
                    alpha_max = alpha = 0.999 * alpha;

                    // Update the objective value at the new trial step
                    x_alpha = x + alpha*dx;
                    f_alpha_max = f_alpha = problem.objective(x_alpha);

                    continue;
                }

                // Check for the Wolfe condition for sufficient decrease in the objective function
                if(f_alpha.val >= f.val + line_search_wolfe*alpha*dot(f.grad, dx))
                {
                    // Decrease the step length
                    alpha *= 0.5;

                    // Update the objective value at the new trial step
                    x_alpha = x + alpha*dx;
                    f_alpha = problem.objective(x_alpha);
                }
            }

            // Check if an adequate step has been found, otherwise use the maximum allowed step length
            if(i >= line_search_max_iterations)
            {
                // Set the step length to its maximum allowed value without causing primal infeasibility
                alpha = std::min(alpha_max, 1.0);

                // Set f(x + alpha*dx) at the maximum allowed step step length
                x_alpha = x + alpha*dx;
                f_alpha = f_alpha_max;
            }

            // Update the primal variables `x` and the dual variables `z`
            x = x_alpha;
            z = -t;

            // Update the result of the objective function at the new iterate `x = x0 + alpha*dx`
            f = f_alpha;
        };

        // Update the error norms
        auto update_errors = [&]()
        {
            // Calculate the current error of the minimization calculation
            error = norm<1>(diag(w)*t);

            // Calculate the feasibility residual
            infeasibility = norminf(A*x - b);
        };

        // Return true if the calculation has converged
        auto converged = [&]()
        {
            Vector tmp = (w.array() > 0).select(z, 0.0);

            if(error < tolerance && min(tmp) > negative_dual_tolerance)
            {
                result.iterations = iter;
                result.succeeded = true;
                return true;
            }
            return false;
        };

        // Perform the minimization of the objective function under linear and bound constraints
        auto minimize = [&]()
        {
            calculate_feasible_point();

            for(; iter < max_iterations; ++iter)
            {
                calculate_descent_direction();

                solve_line_search_minimization_problem();

                update_errors();

                output_state();

                if(converged())
                    break;
            }
        };

        // Return true if there are active primal variables that need to be removed from the boundary
        auto some_active_variables_should_become_inactive = [&]()
        {
            // Check all components in `z` for sufficient negative values
            // and set the corresponding primal variable to an interior value
            bool any_variable_to_become_inactive = false;
            for(unsigned i = 0; i < n; ++i)
            {
                if(z[i] < negative_dual_tolerance)
                {
                    x[i] = l[i] + active_to_inactive;
                    any_variable_to_become_inactive = true;
                }
            }
            return any_variable_to_become_inactive;
        };

        initialize();

        output_header();

        do
        {
            minimize();

        } while(some_active_variables_should_become_inactive());

        // Output the header of the calculation
        outputter.outputHeader();

        // Finish timing the calculation
        result.time = elapsed(begin);

        return result;
    }
};

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

auto OptimumSolverKarpov::solve(const OptimumProblem& problem, OptimumState& state) -> OptimumResult
{
    OptimumOptions options;
    return pimpl->solve(problem, state, options);
}

auto OptimumSolverKarpov::solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    return pimpl->solve(problem, state, options);
}

auto OptimumSolverKarpov::dxdp(const Vector& dgdp, const Vector& dbdp) -> Vector
{
    RuntimeError("Could not calculate the sensitivity of the optimal solution with respect to parameters.",
        "The method OptimumSolverKarpov::dxdp has not been implemented yet.");
    return {};
}

auto OptimumSolverKarpov::clone() const -> OptimumSolverBase*
{
    return new OptimumSolverKarpov(*this);
}

} // namespace Reaktoro
