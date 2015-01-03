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

#include "OptimumSolverIpopt.hpp"

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Common/Outputter.hpp>
#include <Reaktor/Common/TimeUtils.hpp>
#include <Reaktor/Math/MathUtils.hpp>
#include <Reaktor/Optimization/Filter.hpp>
#include <Reaktor/Optimization/KktSolver.hpp>
#include <Reaktor/Optimization/OptimumProblem.hpp>
#include <Reaktor/Optimization/OptimumOptions.hpp>
#include <Reaktor/Optimization/OptimumResult.hpp>
#include <Reaktor/Optimization/Utils.hpp>

namespace Reaktor {

struct OptimumSolverIpopt::Impl
{
    /// The value of the objective function evaluated at the primal solution `x`
    double f;

    /// The gradient of the objective function evaluated at the primal solution `x`
    Vector g;

    /// The value of the equality constraint function evaluated at the primal solution `x`
    Vector h;

    /// The gradient of the equality constraint function evaluated at the primal solution `x`
    Matrix A;

    /// The optimisation filter
    Filter filter;

    /// The components for the the solution of the KKT equations
    Matrix H;
    Matrix invH;
    Vector diagH;
    Vector a, b;
    KktSolver kkt;

    /// The outputter instance
    Outputter outputter;

    /// The auxiliary objective and barrier function results
    double f_trial, phi;
    Vector grad_phi;

    /// The auxiliary constraint function results
    Vector h_trial;

    /// The trial primal iterate
    Vector x_trial;

    /// The Newton steps along x, y and z
    Vector dx, dy, dz;

    /// The trial second-order iterate
    Vector x_soc;

    /// The Newton steps for the trial second-order iterates x and y
    Vector dx_cor, dy_cor;

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

    auto solve(const OptimumProblem& problem, OptimumResult& result, const OptimumOptions& options) -> void;
};

auto OptimumSolverIpopt::Impl::solve(const OptimumProblem& problem, OptimumResult& result, const OptimumOptions& options) -> void
{
    // Start timing the calculation
    Time begin = time();

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
    const unsigned n = problem.numVariables();
    const unsigned m = problem.numConstraints();

    // Initialise the barrier parameter
    double mu = options.ipopt.mu[0];

    // Define some auxiliary references to variables
    auto& x = result.solution.x;
    auto& y = result.solution.y;
    auto& z = result.solution.z;

    // The alpha step sizes used to restric the steps inside the feasible domain
    double alphax, alphaz, alpha;

    // The optimality, feasibility, centrality and total error variables
    double errorf, errorh, errorc, error;

    // The statistics of the calculation
    OptimumStatistics statistics;

    // Ensure the initial guesses for `x` and `y` have adequate dimensions
    if(x.size() != n) x = zeros(n);
    if(y.size() != m) y = zeros(m);
    if(z.size() != n) z = zeros(n);

    // Ensure the initial guess for `x` is inside the feasible domain
    x = max(x, mux*mu*ones(n));

    // Ensure the initial guess for `z` is inside the feasible domain
    z = (z.array() > 0).select(z, mu/x);

    // The transpose representation of matrix `A`
    const auto At = A.transpose();

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
        outputter.addValue(statistics.num_iterations);
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
        g = problem.objectiveGrad(x);
        h = problem.constraint(x);
        A = problem.constraintGrad(x);
    };

    // The function that computes the Newton step based on a regular/dense Hessian scheme
    auto decompose_kkt_regular_hessian = [&]()
    {
        H = problem.objectiveHessian(x);
        H.diagonal() += z/x;
        kkt.decompose(H, A);
    };

    // The function that decomposes the KKT equation based on a diagonal Hessian scheme
    auto decompose_kkt_diagonal_hessian = [&]()
    {
        diagH = problem.objectiveDiagonalHessian(x);
        diagH += z/x;
        kkt.decomposeWithDiagonalH(diagH, A);
    };

    // The function that decomposes the KKT equation based on an inverse Hessian scheme (quasi-Newton approach)
    auto decompose_kkt_inverse_hessian = [&]()
    {
        invH = problem.objectiveInverseHessian(x, g);
        invH = inverseShermanMorrison(invH, z/x);
        kkt.decomposeWithInverseH(invH, A);
    };

    // The function that computes the Newton step
    auto compute_newton_step = [&]()
    {
        // Pre-decompose the KKT equation based on the Hessian scheme
        switch(problem.hessianScheme())
        {
        case HessianScheme::Diagonal: decompose_kkt_diagonal_hessian(); break;
        case HessianScheme::Inverse: decompose_kkt_inverse_hessian(); break;
        default: decompose_kkt_regular_hessian(); break;
        }

        // Compute the right-hand side vectors of the KKT equation
        a.noalias() = -(g - At*y - mu/x);
        b.noalias() = -h;

        // Compute `dx` and `dy` by solving the KKT equation
        kkt.solve(a, b, dx, dy);

        // Compute `dz` with the already computed `dx`
        dz = (mu - z % dx)/x - z;

        statistics.time_linear_system += kkt.statistics().time;
    };

    auto successful_second_order_correction = [&]() -> bool
    {
        h_soc = alpha*h + h_trial;

        double theta_soc_old = theta;

        for(unsigned soc_iter = 0; soc_iter < max_iters_soc; ++soc_iter)
        {
            outputter.outputMessage("...applying the second-order correction step");
            b.noalias() = -h_soc;
            kkt.solve(a, b, dx_cor, dy_cor);
            statistics.time_linear_system += kkt.statistics().time;

            const double alpha_soc = fractionToTheBoundary(x, dx_cor, tau);

            x_soc = x + alpha_soc * dx_cor;

            f_trial = problem.objective(x_soc);
            h_trial = problem.constraint(x_soc);

            // Compute the second-order corrected \theta and \phi measures at the trial iterate
            const double theta_soc = norminf(h_trial);
            const double phi_soc = f_trial - mu * sum(log(x_soc));

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
            const double beta_phi_soc   = phi - gamma_phi * theta_soc;

            // Check if the sufficiency decrease condition is satisfied - the condition (18) of the paper
            const bool sufficient_decrease_soc = lessThan(theta_soc, beta_theta_soc, theta_soc) or lessThan(phi_soc, beta_phi_soc, phi_soc);

            // Check if the Armijo condition is satisfied - the condition (20) of the paper
            const bool armijo_condition_soc = lessThan(phi_soc - phi, eta_phi*alpha_soc*grad_phi_dx, phi);

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
                dx = dx_cor;
                dy = dy_cor;

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
        phi         = f - mu * sum(log(x));
        grad_phi    = g - mu/x;
        grad_phi_dx = dot(grad_phi, dx);
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
        alphax = fractionToTheBoundary(x, dx, tau);
        alphaz = fractionToTheBoundary(z, dz, tau);
        alpha  = alphax;

        // Start the backtracking line-search algorithm
        for(unsigned linesearch_iter = 0; ; ++linesearch_iter)
        {
            // Calculate the trial iterate
            x_trial = x + alpha*dx;

            // Update the objective and constraint states with the trial iterate
            f_trial = problem.objective(x_trial);
            h_trial = problem.constraint(x_trial);

            // Update the barrier objective function with the trial iterate
            const double phi_trial = f_trial - mu * sum(log(x_trial));

            // Compute theta at the trial iterate
            const double theta_trial = norminf(h_trial);

            // Check if the new iterate is acceptable in the filter
            const bool filter_acceptable = filter.acceptable({theta_trial, phi_trial});

            // Calculate some auxiliary variables for checking sufficiency decrease in the measures \theta and \phi
            const double beta_theta = (1 - gamma_theta) * theta;
            const double beta_phi   = phi - gamma_phi * theta;

            // Check if the sufficiency decrease condition is satisfied - the condition (18) of the paper
            sufficient_decrease = lessThan(theta_trial, beta_theta, theta_trial) or lessThan(phi_trial, beta_phi, phi_trial);

            // Check if the switching condition is satisfied - the condition (19) of the paper
            switching_condition = grad_phi_dx < 0 and greaterThan(alpha*pow_phi, delta*pow_theta, alpha*pow_phi);

            // Check if the Armijo condition is satisfied - the condition (20) of the paper
            armijo_condition = lessThan(phi_trial - phi, eta_phi*alpha*grad_phi_dx, phi);

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
        y += alpha * dy;
        z += alphaz * dz;

        // Update the objective and constraint states
        f = f_trial;
        h = h_trial;
        g = problem.objectiveGrad(x);
        A = problem.constraintGrad(x);
    };

    // The function that computes the current error norms
    auto update_errors = [&]()
    {
        // Calculate the optimality, feasibility and centrality errors
        errorf = norminf(g - A.transpose()*y - z);
        errorh = norminf(h);
        errorc = norminf((x % z)/mu - 1);

        // Calculate the maximum error
        error = std::max({errorf, errorh, errorc});
        statistics.error = error;
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
            ++statistics.num_iterations;
        } while(statistics.num_iterations < options.max_iterations);

        // Check if the solution of the subproblem with fixed mu converged
        statistics.converged = statistics.num_iterations < options.max_iterations;
    };

    update_state();
    output_header();

    for(double val : options.ipopt.mu)
    {
        solve_subproblem(val);

        if(not statistics.converged) break;
    }

    result.statistics = statistics;

    // Finish timing the calculation
    result.statistics.time = elapsed(begin);
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

auto OptimumSolverIpopt::solve(const OptimumProblem& problem, OptimumResult& result) -> void
{
    solve(problem, result, {});
}

auto OptimumSolverIpopt::solve(const OptimumProblem& problem, OptimumResult& result, const OptimumOptions& options) -> void
{
    pimpl->solve(problem, result, options);
}

} // namespace Reaktor
