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

#include "AlgorithmIpopt.hpp"

// Armadillo includes
#include <armadillo>

// Reaktor includes
#include <Reaktor/Common/Macros.hpp>
#include <Reaktor/Common/Outputter.hpp>
#include <Reaktor/Optimization/AlgorithmUtils.hpp>
#include <Reaktor/Optimization/SaddlePointUtils.hpp>

namespace Reaktor {

auto ipfeasible(const OptimumProblem& problem, OptimumResult& result, const OptimumOptions& options) -> void
{
    const unsigned n = problem.numVariables();
    const unsigned m = problem.numConstraints();
    result.solution.x  = arma::ones(n);
    result.solution.y  = arma::zeros(m);
    result.solution.zl = arma::zeros(n);
    result.solution.zu = arma::zeros(n);
}

auto ipopt(const OptimumProblem& problem, OptimumResult& result, const OptimumOptions& options) -> void
{
    const auto& objective     = problem.objective();
    const auto& constraint    = problem.constraint();
    const auto& tolerance     = options.tolerance;
    const auto& delta         = options.ipopt.delta;
    const auto& eta_phi       = options.ipopt.eta_phi;
    const auto& gamma_alpha   = options.ipopt.gamma_alpha;
    const auto& gamma_phi     = options.ipopt.gamma_phi;
    const auto& gamma_theta   = options.ipopt.gamma_theta;
    const auto& kappa_soc     = options.ipopt.kappa_soc;
    const auto& max_iters_soc = options.ipopt.max_iters_soc;
    const auto& mu            = options.ipopt.mu;
    const auto& s_phi         = options.ipopt.s_phi;
    const auto& s_theta       = options.ipopt.s_theta;
    const auto& tau           = options.ipopt.tau_min;

    Vector& x = result.solution.x;
    Vector& y = result.solution.y;
    Vector& z = result.solution.zl;

    ObjectiveResult f, f_trial, phi;
    ConstraintResult h, h_trial;

    Vector x_trial;
    Vector dx, dy, dz;

    Vector x_soc;
    Vector dx_cor, dy_cor;

    f = objective(x);
    h = constraint(x);

    const double theta0 = arma::norm(h.func, "inf");
    const double theta_min = 1.0e-4 * std::max(1.0, theta0);
    const double theta_max = 1.0e+4 * std::max(1.0, theta0);

    Filter filter;
    extend(filter, {theta_max, INFINITY});

    SaddlePointProblem saddle_point_problem;
    SaddlePointResult saddle_point_result;

    OptimumStatistics statistics;

    Outputter outputter;

    outputter.setOptions(options.output);

    const unsigned n = problem.numVariables();
    const unsigned m = problem.numConstraints();

    outputter.addEntry("iter");
    outputter.addEntries("x", n);
    outputter.addEntries("y", m);
    outputter.addEntries("z", n);
    outputter.addEntry("f(x)");
    outputter.addEntry("h(x)");
    outputter.addEntry("mu(w)");
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
    outputter.addValue(f.func);
    outputter.addValue(arma::norm(h.func, "inf"));
    outputter.addValue(arma::dot(x, z)/n);
    outputter.addValue("---");
    outputter.addValue("---");
    outputter.addValue("---");
    outputter.addValue("---");
    outputter.addValue("---");
    outputter.addValue("---");
    outputter.addValue("---");
    outputter.outputState();

    do
    {
        phi.func = f.func - mu * std::log(arma::prod(x));
        phi.grad = f.grad - mu/x;

        saddle_point_problem.A = h.grad;
        saddle_point_problem.H = f.hessian + arma::diagmat(z/x);
        saddle_point_problem.f = -(phi.grad - h.grad.t()*y);
        saddle_point_problem.g = -h.func;

        solveNullspace(saddle_point_problem, saddle_point_result);

        dx = saddle_point_result.solution.x;
        dy = saddle_point_result.solution.y;
        dz = -(z % dx + x % z - mu)/x;

        // Calculate some auxiliary variables for calculating alpha_min
        const double theta       = arma::norm(h.func, "inf");
        const double grad_phi_dx = arma::dot(phi.grad, dx);
        const double pow_theta   = std::pow(theta, s_theta);
        const double pow_phi     = std::pow(-grad_phi_dx, s_phi);

        // Calculate the minimum allowed step size alpha_min
        double alpha_min = 0.0;
        if(grad_phi_dx < 0 and theta <= theta_min)
            alpha_min = gamma_alpha * std::min({gamma_theta, -gamma_phi*theta/grad_phi_dx, delta*pow_theta/pow_phi});
        if(grad_phi_dx < 0 and theta > theta_min)
            alpha_min = gamma_alpha * std::min(gamma_theta, -gamma_phi*theta/grad_phi_dx);
        else alpha_min = gamma_alpha * gamma_theta;

        //------------------------------------------------------------------------------
        // Start the backtracking line-search algorithm
        //------------------------------------------------------------------------------
        const double alphax = fractionToTheBoundary(x, dx, tau);
        const double alphaz = fractionToTheBoundary(z, dz, tau);

        double alpha = alphax;

        for(unsigned linesearch_iter = 0; ; ++linesearch_iter)
        {
            // Calculate the trial iterate
            x_trial = x + alpha*dx;

            // Update the objective and constraint states with the trial iterate
            f_trial = objective(x_trial);
            h_trial = constraint(x_trial);

            // Update the barrier objective function with the trial iterate
            const double phi_trial = f_trial.func - mu * std::log(arma::prod(x_trial));

            // Compute theta at the trial iterate
            const double theta_trial = arma::norm(h_trial.func, "inf");

            // Check if the new iterate is acceptable in the filter
            const bool filter_acceptable = acceptable(filter, {theta_trial, phi_trial});

            // Calculate some auxiliary variables for checking sufficiency decrease in the measures \theta and \phi
            const double beta_theta = (1 - gamma_theta) * theta;
            const double beta_phi   = phi.func - gamma_phi * theta;

            // Check if the sufficiency decrease condition is satisfied - the condition (18) of the paper
            const bool sufficient_decrease = lessThan(theta_trial, beta_theta, theta_trial) or lessThan(phi_trial, beta_phi, phi_trial);

            // Check if the switching condition is satisfied - the condition (19) of the paper
            const bool switching_condition = grad_phi_dx < 0 and greaterThan(alpha*pow_phi, delta*pow_theta, alpha*pow_phi);

            // Check if the Armijo condition is satisfied - the condition (20) of the paper
            const bool armijo_condition = lessThan(phi_trial - phi.func, eta_phi*alpha*grad_phi_dx, phi.func);

            const bool theta_condition = lessThan(theta, theta_min, theta);

            // The condition cases of the ipopt algorithm for checking sufficient decrease with respect to the current iterate
            const bool case1_satisfied = ( theta_condition and switching_condition) and armijo_condition;
            const bool case2_satisfied = (!theta_condition or !switching_condition) and sufficient_decrease;

            // Check if the current trial iterate can be accepted
            if(filter_acceptable and (case1_satisfied or case2_satisfied))
            {
                // Extend the filter if the switching condition or the Armijo condition do not hold
                if(switching_condition == false or armijo_condition == false)
                    extend(filter, {beta_theta, beta_phi});

                goto accept_trial_iterate;
            }

            //------------------------------------------------------------------------------
            // Start the second-order correction algorithm
            //------------------------------------------------------------------------------
            if(linesearch_iter == 0 and theta_trial > theta)
            {
                outputter.outputMessage("...applying the second-order correction step");

                Vector hsoc = alpha*h.func + h_trial.func;

                double theta_soc_old = theta;

                for(unsigned soc_iter = 0; soc_iter < max_iters_soc; ++soc_iter)
                {
                    saddle_point_problem.g = -hsoc;
                    solveNullspace(saddle_point_problem, saddle_point_result);
                    dx_cor = saddle_point_result.solution.x;
                    dy_cor = saddle_point_result.solution.y;

                    const double alpha_soc = fractionToTheBoundary(x, dx_cor, tau);

                    x_soc = x + alpha_soc * dx_cor;

                    f_trial = objective(x_soc);
                    h_trial = constraint(x_soc);

                    // Compute the second-order corrected \theta and \phi measures at the trial iterate
                    const double theta_soc = arma::norm(h_trial.func, "inf");
                    const double phi_soc = f_trial.func - mu * std::log(arma::prod(x_soc));

                    // Leave the second-order correction algorithm if \theta is not decreasing sufficiently
                    if(soc_iter > 0 and theta_soc > kappa_soc * theta_soc_old)
                        break;

                    // Leave the second-order correction algorithm if x_soc is not acceptable in the filter
                    if(not acceptable(filter, {theta_soc, phi_soc}))
                        break;

                        // Check if the new iterate is acceptable in the filter
                    const bool filter_acceptable_soc = acceptable(filter, {theta_soc, phi_soc});

                    // Calculate some auxiliary variables for checking sufficiency decrease in the measures \theta and \phi
                    const double beta_theta_soc = (1 - gamma_theta) * theta_soc;
                    const double beta_phi_soc   = phi.func - gamma_phi * theta_soc;

                    // Check if the sufficiency decrease condition is satisfied - the condition (18) of the paper
                    const bool sufficient_decrease_soc = lessThan(theta_soc, beta_theta_soc, theta_soc) or lessThan(phi_soc, beta_phi_soc, phi_soc);

                    // Check if the Armijo condition is satisfied - the condition (20) of the paper
                    const bool armijo_condition_soc = lessThan(phi_soc - phi.func, eta_phi*alpha_soc*grad_phi_dx, phi.func);

                    // The condition cases of the ipopt algorithm for checking sufficient decrease with respect to the current iterate
                    const bool case1_satisfied_soc = (theta <= theta_min and switching_condition) and armijo_condition_soc;
                    const bool case2_satisfied_soc = (theta  > theta_min or !switching_condition) and sufficient_decrease_soc;

                    // Check if the current trial soc iterate can be accepted
                    if(filter_acceptable_soc and (case1_satisfied_soc or case2_satisfied_soc))
                    {
                        // Extend the filter if the switching condition or the Armijo condition do not hold
                        if(switching_condition == false or armijo_condition_soc == false)
                            extend(filter, {beta_theta_soc, beta_phi_soc});

                        // Set the step-size length to the one obtained in the second-order correction algorithm
                        alpha = alpha_soc;

                        // Set the trial iterate to the second-order corrected iterate
                        x_trial = x_soc;

                        // Set the steps dx and dy to their second-order corrected states
                        dx = dx_cor;
                        dy = dy_cor;

                        goto accept_trial_iterate;
                    }

                    // Update the constraint evaluation h^soc
                    hsoc = alpha_soc*hsoc + h_trial.func;

                    // Update the old second-order corrected \theta
                    theta_soc_old = theta_soc;
                }
            }

            // Decrease the length of the size step \alpha
            alpha *= 0.5;

            // Check if the step size is smaller than the minimum
            Assert(alpha > alpha_min, "");
        }

        //------------------------------------------------------------------------------
        // Accept the trial iterate obtained in the backtracking line-search algorithm
        //------------------------------------------------------------------------------
        accept_trial_iterate:

        // Update the iterates x, y, z
        x  = x_trial;
        y +=  alpha * dy;
        z += alphaz * dz;

        // Update the objective and constraint states
        f = f_trial;
        h = h_trial;

        // Calculate the optimality, feasibility and centrality errors
        const double errorf = arma::norm(f.grad - h.grad.t()*y - z, "inf");
        const double errorh = arma::norm(h.func, "inf");
        const double errorc = arma::norm(x%z - mu, "inf");

        // Calculate the maximum error
        statistics.error = std::max({errorf, errorh, errorc});

        ++statistics.num_iterations;

        outputter.addValue(statistics.num_iterations);
        outputter.addValues(x);
        outputter.addValues(y);
        outputter.addValues(z);
        outputter.addValue(f.func);
        outputter.addValue(arma::norm(h.func, "inf"));
        outputter.addValue(arma::dot(x, z)/n);
        outputter.addValue(errorf);
        outputter.addValue(errorh);
        outputter.addValue(errorc);
        outputter.addValue(statistics.error);
        outputter.addValue(alpha);
        outputter.addValue(alphax);
        outputter.addValue(alphaz);
        outputter.outputState();

    } while(statistics.error > tolerance and statistics.num_iterations < options.max_iterations);

    outputter.outputHeader();

    if(statistics.num_iterations < options.max_iterations)
        statistics.converged = true;

    result.statistics = statistics;
}

} // namespace Reaktor
