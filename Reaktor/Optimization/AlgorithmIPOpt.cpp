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

#include "AlgorithmIPOpt.hpp"

// Armadillo includes
#include <armadillo>

// Reaktor includes
#include <Reaktor/Common/Macros.hpp>
#include <Reaktor/Optimization/AlgorithmUtils.hpp>
#include <Reaktor/Optimization/OptimumUtils.hpp>
#include <Reaktor/Optimization/SaddlePointUtils.hpp>

namespace Reaktor {

auto ipfeasible(const OptimumProblem& problem, OptimumResult& result, const OptimumOptions& options) -> void
{
    const unsigned n = problem.numVariables();
    const unsigned m = problem.numConstraints();
    const Vector& l = problem.lowerBounds();
    const Vector& u = problem.upperBounds();
    result.solution.x  = arma::ones(n) + l;
    result.solution.y  = arma::zeros(m);
    result.solution.zl = options.mu/(result.solution.x - l);
    result.solution.zu = options.mu/(u - result.solution.x);
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
//    const auto& kappa_epsilon = options.ipopt.kappa_epsilon;
//    const auto& kappa_mu      = options.ipopt.kappa_mu;
//    const auto& kappa_sigma   = options.ipopt.kappa_sigma;
    const auto& kappa_soc     = options.ipopt.kappa_soc;
    const auto& ksi_phi       = options.ipopt.ksi_phi;
    const auto& max_iters_soc = options.ipopt.max_iters_soc;
    const auto& mu            = options.ipopt.mu;
    const auto& s_phi         = options.ipopt.s_phi;
    const auto& s_theta       = options.ipopt.s_theta;
    const auto& tau_min       = options.ipopt.tau_min;
//    const auto& theta_mu      = options.ipopt.theta_mu;

    Vector& x = result.solution.x;
    Vector& y = result.solution.y;
    Vector& z = result.solution.zl;

    ObjectiveResult f, f_trial, phi, phi_trial;
    ConstraintResult h, h_trial;

    Vector x_trial;
    Vector dx, dy, dz;

    Vector x_soc;
    Vector dx_cor, dy_cor;

    f = objective(x);
    h = constraint(x);

    phi.func = f.func - mu * std::log(arma::prod(x));
    phi.grad = f.grad - mu/x;

    const double theta0 = arma::norm(h.func, "inf");
    const double theta_min = 1.0e-4 * std::max(1.0, theta0);
    const double theta_max = 1.0e+4 * std::max(1.0, theta0);

    Filter filter;
    extend(filter, {theta_max, INFINITY});

    SaddlePointProblem saddle_point_problem;
    SaddlePointResult saddle_point_result;

    do
    {
        // todo maybe there should have a convergence test here to avoid an
        // iteration in case the solution is already provided

        saddle_point_problem.A = h.grad;
        saddle_point_problem.H = f.hessian + arma::diagmat(z/x);
        saddle_point_problem.f = -(phi.grad - h.grad.t()*y);
        saddle_point_problem.g = -h.func;

        solveNullspace(saddle_point_problem, saddle_point_result);

        dx = saddle_point_result.solution.x;
        dy = saddle_point_result.solution.y;
        dz = mu/x - z - (z/x) % dx;

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
        const double tau = std::max(tau_min, 1.0 - mu);
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
            phi_trial.func = f_trial.func - mu * std::log(arma::prod(x_trial));

            // Compute theta at the trial iterate
            const double theta_trial = arma::norm(h_trial.func, "inf");

            // Check if the new iterate is acceptable in the filter
            const bool filter_acceptable = acceptable(filter, {theta_trial, phi_trial.func});

            // Determine some auxiliary Armijo conditions
            const bool armijo_condition1 = theta <= theta_min;
            const bool armijo_condition2 = grad_phi_dx < 0 and alpha*pow_phi > delta*pow_theta;
            const bool armijo_condition3 = lessThan(phi_trial.func - phi.func, eta_phi*alpha*grad_phi_dx, phi.func);

            // Check if the previous Armijo conditions are satisfied
            const bool armijo_satisfied = armijo_condition1 and armijo_condition2 and armijo_condition3;

            // Calculate some auxiliary variables for checking sufficiency decrease in the measures \theta and \phi
            const double beta_theta = (1 - gamma_theta) * theta;
            const double beta_phi   = phi.func - gamma_phi * theta;

            // Check if the current trial iterate sufficiently decreases the measures \theta and \psi
            const bool sufficient_decrease =
                lessThan(theta_trial, beta_theta, theta) or
                lessThan(phi_trial.func, beta_phi, phi.func);

            // Check if the current trial iterate is to be accepted based on the previous determined conditions
            if((filter_acceptable and armijo_satisfied) or (filter_acceptable and sufficient_decrease))
            {
                // Check if the filter needs to be augmented
                if(armijo_condition2 == false or armijo_condition3 == false)
                    extend(filter, {beta_theta, beta_phi});

                goto accept_trial_iterate;
            }

            //------------------------------------------------------------------------------
            // Start the second-order correction algorithm
            //------------------------------------------------------------------------------
            if(linesearch_iter == 0 and theta_trial > theta)
            {
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

                    // Determine some auxiliary Armijo conditions
                    const bool armijo_condition_soc = phi_soc - phi.func - ksi_phi*std::abs(phi.func) <= eta_phi*alpha*grad_phi_dx;

                    // Check if the previous Armijo conditions are satisfied
                    const bool armijo_satisfied_soc = armijo_condition1 and armijo_condition2 and armijo_condition_soc;

                    // Check if the current corrected iterate sufficiently decreases the measures \theta and \psi
                    const bool sufficient_decrease_soc = theta_soc <= beta_theta or phi_soc <= beta_phi;

                    // Check if the current the second-order corrected trial iterate is to be accepted
                    if(armijo_satisfied_soc or sufficient_decrease_soc)
                    {
                        // Set the step-size length to the one obtained in the second-order correction algorithm
                        alpha = alpha_soc;

                        // Set the steps dx and dy to their second-order corrected states
                        dx = dx_cor;
                        dy = dy_cor;

                        // Set the trial iterate to the second-order corrected iterate
                        x_trial = x_soc;

                        // Check if the filter needs to be augmented
                        if(armijo_condition2 == false or armijo_condition_soc == false)
                            extend(filter, {beta_theta, beta_phi});

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
            Assert(alpha > alpha_min, "The step length has achieved a very small number.");
        }

        //------------------------------------------------------------------------------
        // Accept the trial iterate obtained in the backtracking line-search algorithm
        //------------------------------------------------------------------------------
        accept_trial_iterate:

        // Update the iterates x, y, z
        x  = x_trial;
        y += alpha * dy;
        z += alphaz * dz;

        // Update the objective and constraint states
        f = f_trial;
        h = h_trial;

        // Calculate the optimality, feasibility and centrality errors
        const double errorf = arma::norm(f.grad - h.grad.t()*y - z, "inf");
        const double errorh = arma::norm(h.func, "inf");
        const double errorc = arma::norm(x % z, "inf");

        // Calculate the maximum error
        result.statistics.error = std::max({errorf, errorh, errorc});

        ++result.statistics.num_iterations;

    } while(result.statistics.error > tolerance and result.statistics.num_iterations < 100);

    if(result.statistics.num_iterations < 100)
        result.statistics.converged = true;
}

} // namespace Reaktor
