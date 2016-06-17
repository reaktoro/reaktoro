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

#include "OptimumSolverIpNewton.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Outputter.hpp>
#include <Reaktoro/Common/TimeUtils.hpp>
#include <Reaktoro/Math/MathUtils.hpp>
#include <Reaktoro/Optimization/KktSolver.hpp>
#include <Reaktoro/Optimization/OptimumProblem.hpp>
#include <Reaktoro/Optimization/OptimumOptions.hpp>
#include <Reaktoro/Optimization/OptimumResult.hpp>
#include <Reaktoro/Optimization/OptimumState.hpp>
#include <Reaktoro/Optimization/Utils.hpp>

namespace Reaktoro {

struct OptimumSolverIpNewton::Impl
{
    /// The right-hand side of the KKT equations
    KktVector rhs;

    /// The solution of the KKT equations
    KktSolution sol;

    /// The KKT solver
    KktSolver kkt;

    /// The trial iterate x
    Vector xtrial;

    /// The outputter instance
    Outputter outputter;

    /// Solve the optimization problem.
    auto solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
    {
        // Start timing the calculation
        Time begin = time();

        // The result of the calculation
        OptimumResult result;

        // Finish the calculation if the problem has no variable
        if(problem.n == 0)
        {
            state = OptimumState();
            result.succeeded = true;
            result.time = elapsed(begin);
            return result;
        }

        // Initialize the outputter instance
        outputter = Outputter();
        outputter.setOptions(options.output);

        // Set the KKT options
        kkt.setOptions(options.kkt);

        // Define some auxiliary references to variables
        auto& x = state.x;
        auto& y = state.y;
        auto& z = state.z;
        auto& f = state.f;

        // The number of variables and equality constraints
        const auto& A = problem.A;
        const auto& b = problem.b;
        const auto& n = problem.A.cols();
        const auto& m = problem.A.rows();

        // Define auxiliary references to general options
        const auto tol = options.tolerance;
        const auto tolx = options.tolerancex;
        const auto maxiters = options.max_iterations;

        // Define some auxiliary references to IpNewton parameters
        const auto mu = options.ipnewton.mu;
        const auto tau = options.ipnewton.tau;

        // Define some auxiliary references to result variables
        auto& error = result.error;
        auto& iterations = result.iterations;
        auto& succeeded = result.succeeded = false;

        // The regularization parameters delta and gamma
        auto gamma = options.regularization.gamma;
        auto delta = options.regularization.delta;

        // Set gamma and delta to mu in case they are zero
        // This provides even further regularization to the problem,
        // as non-zero gamma and delta prevent unbounded primal and dual
        // variables x and y respectively.
        gamma = gamma ? gamma : mu;
        delta = delta ? delta : mu;

        // Ensure the initial guesses for `x` and `y` have adequate dimensions
        if(x.size() != n) x = zeros(n);
        if(y.size() != m) y = zeros(m);
        if(z.size() != n) z = zeros(n);

        // Ensure the initial guesses for `x` and `z` are inside the feasible domain
        x = (x.array() > 0.0).select(x, mu);
        z = (z.array() > 0.0).select(z, 1.0);

        // The transpose representation of matrix `A`
        const auto At = tr(A);

        // The KKT matrix
        KktMatrix lhs(f.hessian, A, x, z, gamma, delta);

        // The optimality, feasibility, centrality and total error variables
        double errorf, errorh, errorc;

        // The function that outputs the header and initial state of the solution
        auto output_initial_state = [&]()
        {
            if(!options.output.active) return;

            outputter.addEntry("Iteration");
            outputter.addEntries(options.output.xprefix, n, options.output.xnames);
            outputter.addEntries(options.output.yprefix, m, options.output.ynames);
            outputter.addEntries(options.output.zprefix, n, options.output.znames);
            outputter.addEntries("r", n, options.output.xnames);
            outputter.addEntry("f(x)");
            outputter.addEntry("Error");
            outputter.addEntry("Optimality");
            outputter.addEntry("Feasibility");
            outputter.addEntry("Centrality");

            outputter.outputHeader();
            outputter.addValue(iterations);
            outputter.addValues(x);
            outputter.addValues(y);
            outputter.addValues(z);
            outputter.addValues(abs(rhs.rx));
            outputter.addValue(f.val);
            outputter.addValue(error);
            outputter.addValue(errorf);
            outputter.addValue(errorh);
            outputter.addValue(errorc);
            outputter.outputState();
        };

        // The function that outputs the current state of the solution
        auto output_state = [&]()
        {
            if(!options.output.active) return;

            outputter.addValue(iterations);
            outputter.addValues(x);
            outputter.addValues(y);
            outputter.addValues(z);
            outputter.addValues(abs(rhs.rx));
            outputter.addValue(f.val);
            outputter.addValue(error);
            outputter.addValue(errorf);
            outputter.addValue(errorh);
            outputter.addValue(errorc);
            outputter.outputState();
        };

        // Return true if the result of a calculation failed
        auto failed = [&](bool succeeded)
        {
            return !succeeded;
        };

        // The function that computes the current error norms
        auto update_residuals = [&]()
        {
            // Compute the right-hand side vectors of the KKT equation
            rhs.rx.noalias() = -(f.grad - At*y - z + gamma*gamma*ones(n));
            rhs.ry.noalias() = -(A*x + delta*delta*y - b);
            rhs.rz.noalias() = -(x % z - mu);

            // Calculate the optimality, feasibility and centrality errors
            errorf = norminf(rhs.rx);
            errorh = norminf(rhs.ry);
            errorc = norminf(rhs.rz);
            error = std::max({errorf, errorh, errorc});
        };

        // The function that initialize the state of some variables
        auto initialize = [&]()
        {
            // Initialize xtrial
            xtrial.resize(n);

            // Evaluate the objective function
            f = problem.objective(x);

            // Update the residuals of the calculation
            update_residuals();
        };

        // The function that computes the Newton step
        auto compute_newton_step = [&]()
        {
            // Update the decomposition of the KKT matrix with update Hessian matrix
            kkt.decompose(lhs);

            // Compute `dx`, `dy`, `dz` by solving the KKT equation
            kkt.solve(rhs, sol);

            // Update the time spent in linear systems
            result.time_linear_systems += kkt.result().time_solve;
            result.time_linear_systems += kkt.result().time_decompose;

            // Perform emergency Newton step calculation as long as steps contains NaN or INF values
            while(!kkt.result().succeeded)
            {
                // Increase the value of the regularization parameter delta
                delta = std::max(delta * 100, 1e-8);

                // Return false if the calculation did not succeeded
                if(delta > 1e-2) return false;

                // Update the residual of the feasibility conditition
                rhs.ry -= -delta*delta*y;

                // Update the decomposition of the KKT matrix with update Hessian matrix
                kkt.decompose(lhs);

                // Compute `dx`, `dy`, `dz` by solving the KKT equation
                kkt.solve(rhs, sol);

                // Update the time spent in linear systems
                result.time_linear_systems += kkt.result().time_solve;
                result.time_linear_systems += kkt.result().time_decompose;
            }

            // Return true if he calculation succeeded
            return true;
        };

        // The aggressive mode for updating the iterates
        auto update_iterates_aggressive = [&]()
        {
            // Calculate the current trial iterate for x
            for(int i = 0; i < n; ++i)
                xtrial[i] = (x[i] + sol.dx[i] > 0.0) ?
                    x[i] + sol.dx[i] : x[i]*(1.0 - tau);

            // Evaluate the objective function at the trial iterate
            f = problem.objective(xtrial);

            // Initialize the step length factor
            double alpha = fractionToTheBoundary(x, sol.dx, tau);

            // The number of tentatives to find a trial iterate that results in finite objective result
            unsigned tentatives = 0;

            // Repeat until f(xtrial) is finite
            while(!isfinite(f) && ++tentatives < 10)
            {
                // Calculate a new trial iterate using a smaller step length
                xtrial = x + alpha * sol.dx;

                // Evaluate the objective function at the trial iterate
                f = problem.objective(xtrial);

                // Decrease the current step length
                alpha *= 0.5;
            }

            // Return false if xtrial could not be found s.t. f(xtrial) is finite
            if(tentatives == 10)
                return false;

            // Update the iterate x from xtrial
            x = xtrial;

            // Update the z-Lagrange multipliers
            for(int i = 0; i < n; ++i)
                z[i] += (z[i] + sol.dz[i] > 0.0) ?
                    sol.dz[i] : -tau * z[i];

            // Update the y-Lagrange multipliers
            y += sol.dy;

            // Return true as found xtrial results in finite f(xtrial)
            return true;
        };

        // The conservative mode for updating the iterates
        auto update_iterates_convervative = [&]()
        {
            // Initialize the step length factor
            double alphax = fractionToTheBoundary(x, sol.dx, tau);
            double alphaz = fractionToTheBoundary(z, sol.dz, tau);
            double alpha = alphax;

            // The number of tentatives to find a trial iterate that results in finite objective result
            unsigned tentatives = 0;

            // Repeat until a suitable xtrial iterate if found such that f(xtrial) is finite
            for(; tentatives < 10; ++tentatives)
            {
                // Calculate the current trial iterate for x
                xtrial = x + alpha * sol.dx;

                // Evaluate the objective function at the trial iterate
                f = problem.objective(xtrial);

                // Leave the loop if f(xtrial) is finite
                if(isfinite(f))
                    break;

                // Decrease alpha in a hope that a shorter step results f(xtrial) finite
                alpha *= 0.01;
            }

            // Return false if xtrial could not be found s.t. f(xtrial) is finite
            if(tentatives == 10)
                return false;

            // Update the iterate x from xtrial
            x = xtrial;

            // Update the z-Lagrange multipliers
            z += alphaz * sol.dz;

            // Update the y-Lagrange multipliers
            y += sol.dy;

            // Return true as found xtrial results in finite f(xtrial)
            return true;
        };

        // The function that performs an update in the iterates
        auto update_iterates = [&]()
        {
            switch(options.ipnewton.step)
            {
            case Aggressive: return update_iterates_aggressive();
            default: return update_iterates_convervative();
            }
        };

        auto converged = [&]()
        {
            // Check if the calculation should stop based on max variation of x
            if(tolx && max(abs(sol.dx)) < tolx)
                return true;

            // Check if the calculation should stop based on optimality condititions
            return error < tol;
        };

        initialize();
        output_initial_state();

        for(iterations = 1; iterations <= maxiters && !succeeded; ++iterations)
        {
            if(failed(compute_newton_step()))
                break;
            if((succeeded = converged()))
                break;
            if(failed(update_iterates()))
                break;
            update_residuals();
            output_state();
        }

        // Output a final header
        outputter.outputHeader();

        // Finish timing the calculation
        result.time = elapsed(begin);

        return result;
    }

    /// Calculate the sensitivity of the optimal solution with respect to parameters.
    auto dxdp(const Vector& dgdp, const Vector& dbdp) -> Matrix
    {
        // Initialize the right-hand side of the KKT equations
        rhs.rx.noalias() = -dgdp;
        rhs.ry.noalias() =  dbdp;
        rhs.rz.fill(0.0);

        // Solve the KKT equations to get the derivatives
        kkt.solve(rhs, sol);

        // Return the calculated sensitivity vector
        return sol.dx;
    }
};

OptimumSolverIpNewton::OptimumSolverIpNewton()
: pimpl(new Impl())
{}

OptimumSolverIpNewton::OptimumSolverIpNewton(const OptimumSolverIpNewton& other)
: pimpl(new Impl(*other.pimpl))
{}

OptimumSolverIpNewton::~OptimumSolverIpNewton()
{}

auto OptimumSolverIpNewton::operator=(OptimumSolverIpNewton other) -> OptimumSolverIpNewton&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto OptimumSolverIpNewton::solve(const OptimumProblem& problem, OptimumState& state) -> OptimumResult
{
    return pimpl->solve(problem, state, {});
}

auto OptimumSolverIpNewton::solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    return pimpl->solve(problem, state, options);
}

auto OptimumSolverIpNewton::dxdp(const Vector& dgdp, const Vector& dbdp) -> Vector
{
    return pimpl->dxdp(dgdp, dbdp);
}

auto OptimumSolverIpNewton::clone() const -> OptimumSolverBase*
{
    return new OptimumSolverIpNewton(*this);
}

} // namespace Reaktoro
