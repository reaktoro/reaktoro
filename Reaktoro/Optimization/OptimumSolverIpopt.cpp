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

#include "OptimumSolverIpopt.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Outputter.hpp>
#include <Reaktoro/Common/TimeUtils.hpp>
#include <Reaktoro/Math/MathUtils.hpp>
#include <Reaktoro/Optimization/Filter.hpp>
#include <Reaktoro/Optimization/KktSolver.hpp>
#include <Reaktoro/Optimization/OptimumProblem.hpp>
#include <Reaktoro/Optimization/OptimumOptions.hpp>
#include <Reaktoro/Optimization/OptimumResult.hpp>
#include <Reaktoro/Optimization/OptimumState.hpp>
#include <Reaktoro/Optimization/Utils.hpp>

namespace Reaktoro {

struct OptimumSolverIpopt::Impl
{
    /// The optimisation filter
    Filter filter;

    /// The components for the the solution of the KKT equations
    KktVector rhs;
    KktSolution sol;
    KktSolver kkt;

    /// The outputter instance
    Outputter outputter;

    /// The auxiliary objective and barrier function results
    ObjectiveResult f_trial, phi;

    /// The auxiliary constraint function results
    Vector h_trial;

    /// The trial primal iterate
    Vector x_trial;

    /// The trial second-order iterate
    Vector x_soc;

    /// The Newton steps for the trial second-order iterates x and y
    KktSolution sol_cor;

    /// The equality constraint function evaluated at x_soc
    Vector h_soc;

    /// Several other auxiliary algorithmic variables computed along the way
    double theta0;
    double theta_min;
    double theta_max;
    double theta;
    double grad_phi_dx;
    double pow_theta;
    double pow_phi;
    double alpha_min;
    double alphax;
    double alphaz;
    double alpha;
    double errorf;
    double errorh;
    double errorc;
    bool sufficient_decrease;
    bool switching_condition;
    bool armijo_condition;
    bool theta_condition;

    auto solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult;
};

auto OptimumSolverIpopt::Impl::solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    // Start timing the calculation
    Time begin = time();

    // Set the KKT options
    kkt.setOptions(options.kkt);

    // The result of the calculation
    OptimumResult result;

    // Auxiliary references to ipopt parameters
    const auto delta          = options.ipopt.delta;
    const auto eta_phi        = options.ipopt.eta_phi;
    const auto gamma_alpha    = options.ipopt.gamma_alpha;
    const auto gamma_phi      = options.ipopt.gamma_phi;
    const auto gamma_theta    = options.ipopt.gamma_theta;
    const auto kappa_epsilon  = options.ipopt.kappa_epsilon;
    const auto kappa_soc      = options.ipopt.kappa_soc;
    const auto max_iters_soc  = options.ipopt.max_iters_soc;
    const auto mux            = options.ipopt.mux;
    const auto s_phi          = options.ipopt.s_phi;
    const auto s_theta        = options.ipopt.s_theta;
    const auto soc            = options.ipopt.soc;
    const auto tau            = options.ipopt.tau_min;

    // Set the number of variables `n` and equality constraints `m`
    const unsigned n = problem.A.cols();
    const unsigned m = problem.A.rows();

    // Initialize the barrier parameter
    double mu = options.ipopt.mu[0];

    // Define some auxiliary references to variables
    auto& x = state.x;
    auto& y = state.y;
    auto& z = state.z;

    const auto& A = problem.A;
    const auto& b = problem.b;

    ObjectiveResult f;

    Vector h;

    // The alpha step sizes used to restric the steps inside the feasible domain
    double alphax, alphaz, alpha;

    // The optimality, feasibility, centrality and total error variables
    double errorf, errorh, errorc, error;

    // Ensure the initial guesses for `x` and `y` have adequate dimensions
    if(x.size() != n) x = zeros(n);
    if(y.size() != m) y = zeros(m);
    if(z.size() != n) z = zeros(n);

    // Ensure the initial guess for `x` is inside the feasible domain
    x = max(x, mux*mu*ones(n));

    // Ensure the initial guess for `z` is inside the feasible domain
    z = (z.array() > 0).select(z, mu/x);

    // The transpose representation of matrix `A`
    const auto At = tr(problem.A);

    // Set the options of the KKT solver
    kkt.setOptions(options.kkt);

    // The function that outputs the header and initial state of the solution
    auto output_header = [&]()
    {
        if(not options.output.active) return;

        outputter.setOptions(options.output);

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

    auto restart = [&](double new_mu)
    {
        mu = new_mu;

        theta0 = norminf(h);
        theta_min = 1.0e-4 * std::max(1.0, theta0);
        theta_max = 1.0e+4 * std::max(1.0, theta0);

        filter.clear();
        filter.extend({theta_max, infinity()});
    };

    // The function that updates the objective and constraint state
    auto update_state = [&]()
    {
        f = problem.objective(x);
        h = A*x - b;
    };

    // The function that computes the Newton step
    auto compute_newton_step = [&]()
    {
        // Pre-decompose the KKT equation based on the Hessian scheme
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

    auto successful_second_order_correction = [&]() -> bool
    {
        h_soc = alpha*h + h_trial;

        double theta_soc_old = theta;

        for(unsigned soc_iter = 0; soc_iter < max_iters_soc; ++soc_iter)
        {
            outputter.outputMessage("...applying the second-order correction step");
            rhs.ry.noalias() = -h_soc;
            kkt.solve(rhs, sol_cor);
            result.time_linear_systems += kkt.result().time_solve;

            const double alpha_soc = fractionToTheBoundary(x, sol_cor.dx, tau);

            x_soc = x + alpha_soc * sol_cor.dx;

            f_trial = problem.objective(x_soc);
            h_trial = A*x_soc - b;

            // Compute the second-order corrected \theta and \phi measures at the trial iterate
            const double theta_soc = norminf(h_trial);
            const double phi_soc = f_trial.val - mu * sum(log(x_soc));

            // Leave the second-order correction algorithm if \theta is not decreasing sufficiently
            if(soc_iter > 0 and theta_soc > kappa_soc * theta_soc_old)
                break;

            // Leave the second-order correction algorithm if x_soc is not acceptable in the filter
            if(not filter.acceptable({theta_soc, phi_soc}))
                break;

            // Check if the new iterate is acceptable to the filter
            const bool filter_acceptable_soc = filter.acceptable({theta_soc, phi_soc});

            // Calculate some auxiliary variables for checking sufficiency decrease in the measures \theta and \phi
            const double beta_theta_soc = (1 - gamma_theta) * theta_soc;
            const double beta_phi_soc   = phi.val - gamma_phi * theta_soc;

            // Check if the sufficiency decrease condition is satisfied - the condition (18) of the paper
            const bool sufficient_decrease_soc = lessThan(theta_soc, beta_theta_soc, theta_soc) or lessThan(phi_soc, beta_phi_soc, phi_soc);

            // Check if the Armijo condition is satisfied - the condition (20) of the paper
            const bool armijo_condition_soc = lessThan(phi_soc - phi.val, eta_phi*alpha_soc*grad_phi_dx, phi.val);

            // The condition cases of the ipopt algorithm for checking sufficient decrease with respect to the current iterate
            const bool case1_satisfied_soc = (theta <= theta_min and switching_condition) and armijo_condition_soc;
            const bool case2_satisfied_soc = (theta  > theta_min or !switching_condition) and sufficient_decrease_soc;

            // Check if the current trial soc iterate can be accepted
            if(filter_acceptable_soc and (case1_satisfied_soc or case2_satisfied_soc))
            {
                // Extend the filter if the switching condition or the Armijo condition do not hold
                if(switching_condition == false or armijo_condition_soc == false)
                    filter.extend({beta_theta_soc, beta_phi_soc});

                // Set the step-size length to the one obtained in the second-order correction algorithm
                alpha = alpha_soc;

                // Set the trial iterate to the second-order corrected iterate
                x_trial = x_soc;

                // Set the steps dx and dy to their second-order corrected states
                sol.dx = sol_cor.dx;
                sol.dy = sol_cor.dy;

                outputter.outputMessage("...succeeded!\n");

                // The second order correction was succesfully applied
                return true;
            }

            // Update the constraint evaluation h^soc
            h_soc = alpha_soc * h_soc + h_trial;

            // Update the old second-order corrected \theta
            theta_soc_old = theta_soc;
        }

        outputter.outputMessage("...failed!\n");

        // The second order correction could not be succesfully applied
        return false;
    };

    auto successful_backtracking_line_search = [&]() -> bool
    {
        // Calculate some auxiliary variables for calculating alpha_min
        phi.val     = f.val - mu * sum(log(x));
        phi.grad    = f.grad - mu/x;
        grad_phi_dx = dot(phi.grad, sol.dx);
        theta       = norminf(h);
        pow_theta   = std::pow(theta, s_theta);
        pow_phi     = grad_phi_dx < 0 ? std::pow(-grad_phi_dx, s_phi) : 0.0;

        // Calculate the minimum allowed step size alpha_min
        alpha_min = 0.0;
        if(grad_phi_dx < 0 and theta <= theta_min)
            alpha_min = gamma_alpha * std::min({gamma_theta, -gamma_phi*theta/grad_phi_dx, delta*pow_theta/pow_phi});
        if(grad_phi_dx < 0 and theta > theta_min)
            alpha_min = gamma_alpha * std::min(gamma_theta, -gamma_phi*theta/grad_phi_dx);
        else alpha_min = gamma_alpha * gamma_theta;

        // Compute the fraction-to-the-boundary step sizes
        alphax = fractionToTheBoundary(x, sol.dx, tau);
        alphaz = fractionToTheBoundary(z, sol.dz, tau);
        alpha  = alphax;

        // Start the backtracking line-search algorithm
        for(unsigned linesearch_iter = 0; ; ++linesearch_iter)
        {
            // Calculate the trial iterate
            x_trial = x + alpha*sol.dx;

            // Update the objective and constraint states with the trial iterate
            f_trial = problem.objective(x_trial);
            h_trial = A*x_trial - b;

            // Update the barrier objective function with the trial iterate
            const double phi_trial = f_trial.val - mu * sum(log(x_trial));

            // Compute theta at the trial iterate
            const double theta_trial = norminf(h_trial);

            // Check if the new iterate is acceptable in the filter
            const bool filter_acceptable = filter.acceptable({theta_trial, phi_trial});

            // Calculate some auxiliary variables for checking sufficiency decrease in the measures \theta and \phi
            const double beta_theta = (1 - gamma_theta) * theta;
            const double beta_phi   = phi.val - gamma_phi * theta;

            // Check if the sufficiency decrease condition is satisfied - the condition (18) of the paper
            sufficient_decrease = lessThan(theta_trial, beta_theta, theta_trial) or lessThan(phi_trial, beta_phi, phi_trial);

            // Check if the switching condition is satisfied - the condition (19) of the paper
            switching_condition = grad_phi_dx < 0 and greaterThan(alpha*pow_phi, delta*pow_theta, alpha*pow_phi);

            // Check if the Armijo condition is satisfied - the condition (20) of the paper
            armijo_condition = lessThan(phi_trial - phi.val, eta_phi*alpha*grad_phi_dx, phi.val);

            // Check if the theta condition is satisfied
            theta_condition = lessThan(theta, theta_min, theta);

            // The condition cases of the ipopt algorithm for checking sufficient decrease with respect to the current iterate
            const bool case1_satisfied = ( theta_condition and switching_condition) and armijo_condition;
            const bool case2_satisfied = (!theta_condition or !switching_condition) and sufficient_decrease;

            // Check if the current trial iterate can be accepted
            if(filter_acceptable and (case1_satisfied or case2_satisfied))
            {
                // Extend the filter if the switching condition or the Armijo condition do not hold
                if(switching_condition == false or armijo_condition == false)
                    filter.extend({beta_theta, beta_phi});
                return true;
            }

            // Check if the second order corrections can be applied and if it succeeds
            if(linesearch_iter == 0 and theta_trial > theta and soc)
                if(successful_second_order_correction())
                    return true;

            // Decrease the length of the size step alpha
            alpha *= 0.5;

            // Check if the step size is smaller than the minimum
            if(alpha < alpha_min)
                return false;
        }
    };

    auto accept_trial_iterate = [&]()
    {
        // Update the next iterate
        x  = x_trial;
        y += alpha * sol.dy;
        z += alphaz * sol.dz;

        // Update the objective and constraint states
        f = f_trial;
        h = h_trial;
    };

    // The function that computes the current error norms
    auto update_errors = [&]()
    {
        // Calculate the optimality, feasibility and centrality errors
        errorf = norminf(f.grad - At*y - z);
        errorh = norminf(h);
        errorc = norminf(x % z - mu);

        // Calculate the maximum error
        error = std::max({errorf, errorh, errorc});
        result.error = error;
    };

    auto converged = [&]() -> bool
    {
        // Calculte the tolerance for the convergence checking
        const double tol = std::max(options.tolerance, kappa_epsilon*mu);
        return error < tol;
    };

    auto solve_subproblem = [&](double mu)
    {
        restart(mu);

        do {
            update_errors();
            output_state();
            if(converged()) break;
            compute_newton_step();
            successful_backtracking_line_search();
            accept_trial_iterate();
            ++result.iterations;
        } while(result.iterations < options.max_iterations);

        // Check if the solution of the subproblem with fixed mu converged
        result.succeeded = result.iterations < options.max_iterations;
    };

    update_state();
    output_header();

    for(double val : options.ipopt.mu)
    {
        solve_subproblem(val);

        if(not result.succeeded) break;
    }

    // Finish timing the calculation
    result.time = elapsed(begin);

    return result;
}

OptimumSolverIpopt::OptimumSolverIpopt()
: pimpl(new Impl())
{}

OptimumSolverIpopt::OptimumSolverIpopt(const OptimumSolverIpopt& other)
: pimpl(new Impl(*other.pimpl))
{}

OptimumSolverIpopt::~OptimumSolverIpopt()
{}

auto OptimumSolverIpopt::operator=(OptimumSolverIpopt other) -> OptimumSolverIpopt&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto OptimumSolverIpopt::solve(const OptimumProblem& problem, OptimumState& state) -> OptimumResult
{
    return solve(problem, state, {});
}

auto OptimumSolverIpopt::solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    return pimpl->solve(problem, state, options);
}

} // namespace Reaktoro
