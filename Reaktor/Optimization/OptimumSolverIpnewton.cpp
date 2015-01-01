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

#include "OptimumSolverIpnewton.hpp"

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Common/Macros.hpp>
#include <Reaktor/Common/Outputter.hpp>
#include <Reaktor/Common/TimeUtils.hpp>
#include <Reaktor/Math/MathUtils.hpp>
#include <Reaktor/Optimization/KktSolver.hpp>
#include <Reaktor/Optimization/OptimumProblem.hpp>
#include <Reaktor/Optimization/OptimumOptions.hpp>
#include <Reaktor/Optimization/OptimumResult.hpp>
#include <Reaktor/Optimization/Utils.hpp>

namespace Reaktor {

struct OptimumSolverIpnewton::Impl
{
    const OptimumProblem* problem;

    OptimumResult* result;

    const OptimumOptions* options;

    double f;
    Vector g;
    Vector h;

    Vector dx, dy, dz;

    KktMatrix lhs;
    KktVector rhs;
    KktResult kkt_result;
    KktProblem kkt_problem;

    Outputter outputter;

    Impl()
    {}

    auto solve(const OptimumProblem& problem, OptimumResult& result, const OptimumOptions& options) -> void
    {
        Time begin = time();

        // Define some auxiliary references to variables
        auto& x = result.solution.x;
        auto& y = result.solution.y;
        auto& z = result.solution.z;

        // Define some auxiliary references to parameters
        const auto& n         = problem.numVariables();
        const auto& m         = problem.numConstraints();
        const auto& tolerance = options.tolerance;
        const auto& mu        = options.ipnewton.mu;
        const auto& mux       = options.ipnewton.mux;
        const auto& tau       = options.ipnewton.tau;

        // Ensure the initial guesses for `x` and `y` have adequate dimensions
        if(x.size() != n) x = zeros(n);
        if(y.size() != m) y = zeros(m);

        // Ensure the initial guess for `x` is inside the feasible domain
        x = max(x, mux*mu*ones(n));

        // Ensure the initial guess for `z` is on the central line
        z = mu/x;

        // The reference to the transpose representation of matrix `A`
        auto& A = lhs.A;
        const auto At = A.transpose();

        // The alpha step sizes used to restric the steps inside the feasible domain
        double alphax, alphaz, alpha;

        // The optimality, feasibility, centrality and total error variables
        double errorf, errorh, errorc, error;

        // The statistics of the calculation
        OptimumStatistics statistics;

        // The function that outputs the header and initial state of the solution
        auto output_header = [&]()
        {
            if(not options.output.active) return;

            outputter.setOptions(options.output);

            outputter.addEntry("iter");
            outputter.addEntries("x", n);
            outputter.addEntries("y", m);
            outputter.addEntries("z", n);
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
            outputter.addValue(result.statistics.num_iterations);
            outputter.addValues(x);
            outputter.addValues(y);
            outputter.addValues(z);
            outputter.addValue(f);
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

            outputter.addValue(result.statistics.num_iterations);
            outputter.addValues(x);
            outputter.addValues(y);
            outputter.addValues(z);
            outputter.addValue(f);
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
            g = problem.objectiveGrad(x);
            h = problem.constraint(x);
            A = problem.constraintGrad(x);
        };

        // The function that computes the Newton step based on a regular/dense Hessian scheme
        auto compute_newton_step_regular_hessian = [&]()
        {
            lhs.H = problem.objectiveHessian(x);

            lhs.H.diagonal() += z/x;
            rhs.f = -(g - At*y - mu/x);
            rhs.g = -h;

            kkt_problem.decompose(lhs);
            kkt_problem.solve(lhs, rhs, kkt_result);
        };

        // The function that computes the Newton step based on a diagonal Hessian scheme
        auto compute_newton_step_diagonal_hessian = [&]()
        {
            lhs.diagH = problem.objectiveDiagonalHessian(x);

            lhs.diagH += z/x;
            rhs.f = -(g - At*y - mu/x);
            rhs.g = -h;

            kkt_problem.decomposeWithDiagonalH(lhs);
            kkt_problem.solveWithDiagonalH(lhs, rhs, kkt_result);
        };

        // The function that computes the Newton step based on a inverse Hessian scheme (quasi-Newton approach)
        auto compute_newton_step_inverse_hessian = [&]()
        {
            lhs.invH = problem.objectiveInverseHessian(x, g);

            // Compute the inverse of the matrix `H + inv(X)Z` with known inv(H)
            lhs.invH = inverseShermanMorrison(lhs.invH, z/x);
            rhs.f = -(g - At*y - mu/x);
            rhs.g = -h;

            kkt_problem.decomposeWithInverseH(lhs);
            kkt_problem.solveWithInverseH(lhs, rhs, kkt_result);
        };

        // The function that computes the Newton step based on the current settings of the Hessian scheme
        auto compute_newton_step = [&]()
        {
            switch(problem.hessianScheme())
            {
            case HessianScheme::Diagonal: compute_newton_step_diagonal_hessian(); break;
            case HessianScheme::Inverse: compute_newton_step_inverse_hessian(); break;
            default: compute_newton_step_regular_hessian(); break;
            }

            dx = kkt_result.solution.x;
            dy = kkt_result.solution.y;
            dz = (mu - z % dx)/x - z;

            statistics.time_linear_system += kkt_result.statistics.time;
        };

        // The function that performs an update in the iterates
        auto update_iterates = [&]()
        {
            alphax = fractionToTheBoundary(x, dx, tau);
            alphaz = fractionToTheBoundary(z, dz, tau);
            alpha  = std::min(alphax, alphaz);

            if(options.ipnewton.uniform_newton_step)
            {
                x += alpha * dx;
                y += alpha * dy;
                z += alpha * dz;
            }
            else
            {
                x += alpha * dx;
                y += dy;
                z += alphaz * dz;
            }
        };

        // The function that computes the current error norms
        auto update_errors = [&]()
        {
            // Calculate the optimality, feasibility and centrality errors
            errorf = norminf(g - At*y - z);
            errorh = norminf(h);
            errorc = norminf(x%z - mu);

            // Calculate the maximum error
            error = std::max({errorf, errorh, errorc});
            statistics.error = error;
        };

        update_state();
        output_header();

        do
        {
            ++statistics.num_iterations;
            compute_newton_step();
            update_iterates();
            update_state();
            update_errors();
            output_state();
        } while(error > tolerance and statistics.num_iterations < options.max_iterations);

        outputter.outputHeader();

        if(result.statistics.num_iterations < options.max_iterations)
            statistics.converged = true;

        result.statistics = statistics;
        result.statistics.time = elapsed(begin);
    }
};

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

auto OptimumSolverIpnewton::solve(const OptimumProblem& problem, OptimumResult& result) -> void
{
    pimpl->solve(problem, result, {});
}

auto OptimumSolverIpnewton::solve(const OptimumProblem& problem, OptimumResult& result, const OptimumOptions& options) -> void
{
    pimpl->solve(problem, result, options);
}


} // namespace Reaktor
