// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "OptimumSolverIpBounds.hpp"

// Eigen includes
#include <Reaktoro/deps/eigen3/Eigen/LU>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Outputter.hpp>
#include <Reaktoro/Common/TimeUtils.hpp>
#include <Reaktoro/Math/MathUtils.hpp>
#include <Reaktoro/Optimization/KktSolver.hpp>
#include <Reaktoro/Optimization/OptimumOptions.hpp>
#include <Reaktoro/Optimization/OptimumProblem.hpp>
#include <Reaktoro/Optimization/OptimumResult.hpp>
#include <Reaktoro/Optimization/OptimumState.hpp>
#include <Reaktoro/Optimization/Utils.hpp>

namespace Reaktoro {

struct OptimumSolverIpBounds::Impl
{
    /// The right-hand side vector of the KKT equations
    Vector rhs;

    /// The left-hand side matrix of the KKT equations
    Matrix lhs;

    /// The slack variables `s` and its inverse
    Vector s, inv_s;

    /// The Newton steps `dx`, `dz`, and `ds`
    Vector dx, dz, ds;

    /// The vectors containing the residual of the optimality and feasibility equations
    Vector res_o, res_f;

    /// The trial iterate x
    Vector xtrial;

    /// The outputter instance
    Outputter outputter;

    /// Solve the optimization problem.
    auto solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
    {
        //        // Start timing the calculation
        //        Time begin = time();
        //
        //        // Initialize the outputter instance
        //        outputter = Outputter();
        //        outputter.setOptions(options.output);
        //
        // The result of the calculation
        OptimumResult result;
        //
        //        // The number of primal variables (n) and inequality constraints (m)
        //        const auto n = problem.Ai.cols();
        //        const auto m = problem.Ai.rows();
        //
        //        // Define some auxiliary references to variables
        //        auto& x = state.x;
        //        auto& z = state.z;
        //        auto& f = state.f;
        //
        //        // Define auxiliary references to problem data
        //        const auto& A  = problem.Ai;
        //        const auto& b  = problem.bi;
        //        const auto& At = tr(A);
        //
        //        // Define auxiliary references to general options
        //        const auto tol = options.tolerance;
        //        const auto tolx = options.tolerancex;
        //        const auto maxiters = options.max_iterations;
        //
        //        // Define some auxiliary references to IpBounds parameters
        //        const auto mu = options.ipnewton.mu;
        //        const auto tau = options.ipnewton.tau;
        //
        //        // Define some auxiliary references to result variables
        //        auto& error = result.error;
        //        auto& iterations = result.iterations;
        //        auto& succeeded = result.succeeded = false;
        //
        //        // Ensure the initial guesses for x and z have adequate dimensions
        //        if(x.size() != n) x = zeros(n);
        //        if(z.size() != m) z = zeros(m);
        //
        //        // Initialize the slack variables `s`
        //        s = A*x - b;
        //
        //        // Ensure the initial guesses for `s` and `z` have positive components
        //        z = (z.array() > 0.0).select(z, 1.0);
        //        s = (s.array() > 0.0).select(s, 1.0);
        //
        //        // The optimality, feasibility, and centrality error variables
        //        double error_o, error_f, error_c;
        //
        //        // The function that outputs the header and initial state of the solution
        //        auto output_initial_state = [&]()
        //        {
        //            if(!options.output.active) return;
        //
        //            outputter.addEntry("Iteration");
        //            outputter.addEntries(options.output.xprefix, n, options.output.xnames);
        //            outputter.addEntries(options.output.zprefix, m, options.output.znames);
        //            outputter.addEntries("r", n, options.output.xnames);
        //            outputter.addEntry("f(x)");
        //            outputter.addEntry("Error");
        //            outputter.addEntry("Optimality");
        //            outputter.addEntry("Feasibility");
        //            outputter.addEntry("Centrality");
        //
        //            outputter.outputHeader();
        //            outputter.addValue(iterations);
        //            outputter.addValues(x);
        //            outputter.addValues(z);
        //            outputter.addValues(abs(rhs));
        //            outputter.addValue(f.val);
        //            outputter.addValue(error);
        //            outputter.addValue(error_o);
        //            outputter.addValue(error_f);
        //            outputter.addValue(error_c);
        //            outputter.outputState();
        //        };
        //
        //        // The function that outputs the current state of the solution
        //        auto output_state = [&]()
        //        {
        //            if(!options.output.active) return;
        //
        //            outputter.addValue(iterations);
        //            outputter.addValues(x);
        //            outputter.addValues(z);
        //            outputter.addValues(abs(rhs));
        //            outputter.addValue(f.val);
        //            outputter.addValue(error);
        //            outputter.addValue(error_o);
        //            outputter.addValue(error_f);
        //            outputter.addValue(error_c);
        //            outputter.outputState();
        //        };
        //
        //        // Return true if the result of a calculation failed
        //        auto failed = [&](bool succeeded)
        //        {
        //            return !succeeded;
        //        };
        //
        //        // The function that computes the current error norms
        //        auto update_residuals = [&]()
        //        {
        //            // Update the inverse of the slack variables `s`
        //            inv_s.noalias() = 1.0/s;
        //
        //            // Update the residual vectors of the optimality and feasibility equations
        //            res_o = f.grad - At*z;
        //            res_f = A*x - s - b;
        //
        //            // Compute the right-hand side vector of the KKT equation
        //            rhs.noalias() = -f.grad + mu*At*inv_s - At*((z/s)%res_f);
        //
        //            // Calculate the optimality, feasibility and centrality errors
        //            error_o = norminf(res_o);
        //            error_f = norminf(res_f);
        //            error_c = norminf(s % z - mu);
        //            error = std::max({error_o, error_f, error_c});
        //        };
        //
        //        // The function that initialize the state of some variables
        //        auto initialize = [&]()
        //        {
        //            // Initialize xtrial
        //            xtrial.resize(n);
        //
        //            // Evaluate the objective function
        //            f = problem.objective(x);
        //
        //            // Update the residuals of the calculation
        //            update_residuals();
        //        };
        //
        //        // The function that computes the Newton step
        //        auto compute_newton_step = [&]()
        //        {
        //            // Assemble the KKT matrix
        //            lhs = f.hessian.dense; lhs += At * diag(z/s) * A;
        //
        //            // Compute the steps dx, dz, ds
        //            dx = -lhs.lu().solve(rhs);
        //            ds = A*(x + dx) - s - b;
        //            dz = mu*inv_s - z - z % ds % inv_s;
        //        };
        //
        //        // The function that performs an update in the iterates
        //        auto update_iterates = [&]()
        //        {
        //            // Initialize the step length factor for Newton step dx with the largest possible value
        //            double alphax = fractionToTheBoundary(x, dx, A, b, tau);
        //
        //            // Initialize the step length factor for the trial calculation of dx below
        //            double alpha = 1.0;
        //
        //            // The number of tentatives to find a trial iterate that results in finite objective result
        //            unsigned tentatives = 0;
        //
        //            // Repeat until a suitable xtrial iterate if found such that f(xtrial) is finite
        //            for(; tentatives < 6; ++tentatives)
        //            {
        //                // Calculate the current trial iterate for x
        //                xtrial = x + alpha*alphax*dx;
        //
        //                // Evaluate the objective function at the trial iterate
        //                f = problem.objective(xtrial);
        //
        //                // Leave the loop if f(xtrial) is finite
        //                if(isfinite(f))
        //                    break;
        //
        //                // Decrease alpha in a hope that a shorter step results f(xtrial) finite
        //                alpha *= 0.1;
        //            }
        //
        //            // Return false if xtrial could not be found s.t. f(xtrial) is finite
        //            if(tentatives == 6)
        //                return false;
        //
        //            // Update the iterate x from xtrial
        //            x = xtrial;
        //
        //            // Update the slack variables s
        //            for(int i = 0; i < m; ++i)
        //                s[i] += (s[i] + alpha*ds[i] > 0.0) ?
        //                    alpha*ds[i] : -alpha*tau * s[i];
        //
        //            // Update the Lagrange multipliers z
        //            for(int i = 0; i < m; ++i)
        //                z[i] += (z[i] + alpha*dz[i] > 0.0) ?
        //                    alpha*dz[i] : -alpha*tau * z[i];
        //
        //            // Return true as found xtrial results in finite f(xtrial)
        //            return true;
        //        };
        //
        //        auto converged = [&]()
        //        {
        //            // Check if the calculation should stop based on max variation of x
        //            if(tolx && max(abs(dx)) < tolx)
        //                return true;
        //
        //            // Check if the calculation should stop based on optimality condititions
        //            return error < tol;
        //        };
        //
        //        initialize();
        //        output_initial_state();
        //
        //        for(iterations = 1; iterations <= maxiters && !succeeded; ++iterations)
        //        {
        //            compute_newton_step();
        //            if(failed(update_iterates()))
        //                break;
        //            update_residuals();
        //            output_state();
        //            succeeded = converged();
        //        }
        //
        //        // Output a final header
        //        outputter.outputHeader();
        //
        //        // Finish timing the calculation
        //        result.time = elapsed(begin);
        //
        return result;
    }

    /// Calculate the sensitivity of the optimal solution with respect to parameters.
    auto dxdp(VectorConstRef dgdp, VectorConstRef dbdp) -> Matrix
    {
        RuntimeError("Could not calculate the sensitivity of the optimal solution with respect to parameters.",
                     "The method OptimumSolverIpBounds::dxdp has not been implemented yet.");
        return {};
    }
};

OptimumSolverIpBounds::OptimumSolverIpBounds()
    : pimpl(new Impl())
{}

OptimumSolverIpBounds::OptimumSolverIpBounds(const OptimumSolverIpBounds& other)
    : pimpl(new Impl(*other.pimpl))
{}

OptimumSolverIpBounds::~OptimumSolverIpBounds()
{}

auto OptimumSolverIpBounds::operator=(OptimumSolverIpBounds other) -> OptimumSolverIpBounds&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto OptimumSolverIpBounds::solve(const OptimumProblem& problem, OptimumState& state) -> OptimumResult
{
    return pimpl->solve(problem, state, {});
}

auto OptimumSolverIpBounds::solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    return pimpl->solve(problem, state, options);
}

auto OptimumSolverIpBounds::dxdp(VectorConstRef dgdp, VectorConstRef dbdp) -> Vector
{
    return pimpl->dxdp(dgdp, dbdp);
}

auto OptimumSolverIpBounds::clone() const -> OptimumSolverBase*
{
    return new OptimumSolverIpBounds(*this);
}

} // namespace Reaktoro
