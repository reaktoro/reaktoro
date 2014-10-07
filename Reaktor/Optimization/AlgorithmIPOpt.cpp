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

auto ipopt(const OptimumProblem& problem, OptimumResult& result, const OptimumOptions& options) -> void
{
    const unsigned n = numVariables(problem);
    const unsigned m = numConstraints(problem);

    const double tolerance = options.tolerance;
    const double mu        = options.mu;

    const double delta           = 1.0;
    const double eta_phi         = 1.0e-4;
    const double gamma_alpha     = 0.05;
    const double gamma_phi       = 1.0e-5;
    const double gamma_theta     = 1.0e-5;
    const double kappa_epsilon   = 10.0;
    const double kappa_mu        = 0.2;
    const double kappa_sigma     = 100;
    const double kappa_soc       = 0.99;
    const double ksi_phi         = 1.0e-15;
    const double s_phi           = 2.3;
    const double s_theta         = 1.1;
    const double tau             = 0.995;
    const double theta_mu        = 1.5;
    const unsigned max_iters_soc = 4;

    Vector& x = result.solution.x;
    Vector& y = result.solution.y;
    Vector& z = result.solution.zl;

    ObjectiveResult f, f_trial;

    Vector h, h_trial;

    Vector x_trial;
    Vector dx, dy, dz;

    Vector x_soc;
    Vector dx_cor, dy_cor;

    // todo perhaps this should be left to the user
    if(x.size() != n) x = arma::ones(n);
    if(y.size() != m) y = arma::zeros(m);
    if(z.size() != n) z = mu/x;

    f = objective(problem, x);
    h = constraint(problem, x);

    const double theta0 = arma::norm(h, "inf");
    const double theta_min = 1.0e-4 * std::max(1.0, theta0);
    const double theta_max = 1.0e+4 * std::max(1.0, theta0);

    Filter filter;
    extend({theta_max, INFINITY}, filter);

    SaddlePointProblem saddle_point_problem;
    SaddlePointResult saddle_point_result;

    const Matrix& A = problem.A;

    do
    {
        const double theta = arma::norm(h, "inf");
        const double phi = f.func - mu * std::log(arma::prod(x));

        Vector invXZ = z/x;
        Vector uXe = mu/x;
        Vector grad_phi = f.grad - uXe;

        saddle_point_problem.A = A;
        saddle_point_problem.H = f.hessian + arma::diagmat(invXZ);
        saddle_point_problem.f = -grad_phi + A.t()*y;
        saddle_point_problem.g = -h;

        solveNullspace(saddle_point_problem, saddle_point_result);

        dx = saddle_point_result.solution.x;
        dy = saddle_point_result.solution.y;
        dz = uXe - z - invXZ % dx;

        const double alphax = std::min(1.0, tau * largestStep(x, dx));
        const double alphaz = std::min(1.0, tau * largestStep(z, dz));

        // Calculate some auxiliary variables for calculating \alpha_{min}
        const double grad_phi_dx = arma::dot(grad_phi, dx);
        const double pow_theta = std::pow(theta, s_theta);
        const double pow_phi = std::pow(-grad_phi_dx, s_phi);

        // Calculate the minimum allowed step size \alpha_{min}
        double alpha_min;

        if(grad_phi_dx < 0 and theta <= theta_min)
            alpha_min = gamma_alpha * std::min({gamma_theta, -gamma_phi*theta/grad_phi_dx, delta*pow_theta/pow_phi});

        if(grad_phi_dx < 0 and theta > theta_min)
            alpha_min = gamma_alpha * std::min(gamma_theta, -gamma_phi*theta/grad_phi_dx);

        else alpha_min = gamma_alpha * gamma_theta;


        alpha_min = mu*mu;

        //------------------------------------------------------------------------------
        // Start the backtracking line-search algorithm
        //------------------------------------------------------------------------------
        double alpha = alphax;

        for(unsigned linesearch_iter = 0; ; ++linesearch_iter)
        {
            // Calculate the trial iterate x(trial)
            x_trial = x + alpha * dx;

            // Update the objective and constraint states with the trial iterates
            f_trial = objective(problem, x_trial);
            h_trial = constraint(problem, x_trial);

            // Compute the \theta and \phi measures at the trial iterate
            const double theta_trial = arma::norm(h_trial, "inf");
            const double phi_trial = f_trial.func - mu * std::log(arma::prod(x_trial));

            // Check if the new iterate is acceptable in the filter
            const bool filter_acceptable = acceptable({theta_trial, phi_trial}, filter);

            // Determine some auxiliary Armijo conditions
            const bool armijo_condition1 = theta_trial <= theta_min;
            const bool armijo_condition2 = grad_phi_dx < 0 and alpha*pow_phi > delta*pow_theta;
            const bool armijo_condition3 = lessThan(phi_trial - phi, eta_phi*alpha*grad_phi_dx, phi);

            // Check if the previous Armijo conditions are satisfied
            const bool armijo_satisfied = armijo_condition1 and armijo_condition2 and armijo_condition3;

            // Calculate some auxiliary variables for checking sufficiency decrease in the measures \theta and \phi
            const double beta_theta = (1 - gamma_theta) * theta;
            const double beta_phi = phi - gamma_phi * theta;

            // Check if the current trial iterate sufficiently decreases the measures \theta and \psi
            const bool sufficient_decrease =
                lessThan(theta_trial, beta_theta, theta) or
                lessThan(phi_trial, beta_phi, phi);

            // Check if the current trial iterate is to be accepted based on the previous determined conditions
            if(filter_acceptable and armijo_satisfied or sufficient_decrease)
            {
                // Check if the filter needs to be augmented
                if(armijo_condition2 == false or armijo_condition3 == false)
                    extend({beta_theta, beta_phi}, filter);

                goto accept_trial_iterate;
            }

            //------------------------------------------------------------------------------
            // Start the second-order correction algorithm
            //------------------------------------------------------------------------------
            if(linesearch_iter == 0 and theta_trial > theta)
            {
                Vector hsoc = alpha*h + h_trial;

                double theta_soc_old = theta;

                for(unsigned soc_iter = 0; soc_iter < max_iters_soc; ++soc_iter)
                {
                    saddle_point_problem.g = -hsoc;
                    solveNullspace(saddle_point_problem, saddle_point_result);
                    dx_cor = saddle_point_result.solution.x;
                    dy_cor = saddle_point_result.solution.y;

                    const double alpha_soc = std::min(1.0, tau * largestStep(x, dx_cor));

                    x_soc = x + alpha_soc * dx_cor;

                    f_trial = objective(problem, x_soc);
                    h_trial = constraint(problem, x_soc);

                    // Compute the second-order corrected \theta and \phi measures at the trial iterate
                    const double theta_soc = arma::norm(h_trial, "inf");
                    const double phi_soc = f_trial.func - mu * std::log(arma::prod(x_soc));

                    // Leave the second-order correction algorithm if \theta is not decreasing sufficiently
                    if(soc_iter > 0 and theta_soc > kappa_soc * theta_soc_old)
                        break;

                    // Leave the second-order correction algorithm if x_soc is not acceptable in the filter
                    if(not acceptable({theta_soc, phi_soc}, filter))
                        break;

                    // Determine some auxiliary Armijo conditions
                    const bool armijo_condition_soc = phi_soc - phi - ksi_phi*std::abs(phi) <= eta_phi*alpha*grad_phi_dx;

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
                            extend({beta_theta, beta_phi}, filter);

                        goto accept_trial_iterate;
                    }

                    // Update the constraint evaluation h^soc
                    hsoc = alpha_soc*hsoc + h_trial;

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
        y +=  alpha * dy;
        z += alphaz * dz;

        // Update the objective and constraint states
        f = f_trial;
        h = h_trial;

        // Calculate the optimality, feasibility and centrality errors
        const double errorf = arma::norm(f.grad - A.t()*y - z, "inf");
        const double errorh = arma::norm(h, "inf");
        const double errorc = arma::norm(x % z, "inf");

        // Calculate the maximum error
        result.statistics.error = std::max({errorf, errorh, errorc});

        ++result.statistics.num_iterations;

    } while(result.statistics.error > tolerance and result.statistics.num_iterations < 100);

    if(result.statistics.num_iterations < 100)
        result.statistics.converged = true;
}

} // namespace Reaktor
