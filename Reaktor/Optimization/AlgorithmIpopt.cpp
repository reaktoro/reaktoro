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
#include <Reaktor/Optimization/FilterUtils.hpp>
#include <Reaktor/Optimization/OptimumOptions.hpp>
#include <Reaktor/Optimization/OptimumProblem.hpp>
#include <Reaktor/Optimization/OptimumResult.hpp>
#include <Reaktor/Optimization/SaddlePointUtils.hpp>

namespace Reaktor {
namespace {

struct IpoptSolver
{
    /// The optimization problem instance
    const OptimumProblem& problem;

    /// The optimization result instance
    OptimumResult& result;

    /// The optimization options instance
    const OptimumOptions& options;

    /// The auxiliary references to the primal and dual solutions
    Vector &x, &y, &zl, &zu;

    /// The auxiliary reference to the statistics of the calculation
    OptimumStatistics& statistics;

    /// The auxiliary references to the lower and upper bounds
    const Vector &lower, &upper;

    const SaddlePointOptions& saddle_point_options;

    ObjectiveResult f, f_trial, phi;
    ConstraintResult h, h_trial;

    Vector x_trial;
    Vector dx, dy, dz;
    Vector x_soc;
    Vector dx_cor, dy_cor;
    Vector hsoc;
    Vector dl, du, d;
    unsigned n, m;
    double mu;
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
    Filter filter;
    SaddlePointProblem saddle_point_problem;
    SaddlePointResult saddle_point_result;
    Outputter outputter;

    IpoptSolver(const OptimumProblem& problem, OptimumResult& result, const OptimumOptions& options)
    : problem(problem), result(result), options(options),
      x(result.solution.x), y(result.solution.y),
      zl(result.solution.zl), zu(result.solution.zu), statistics(result.statistics),
      lower(problem.lowerBounds()), upper(problem.upperBounds()),
      saddle_point_options(options.ipopt.saddle_point),
      n(problem.numVariables()), m(problem.numConstraints())
    {
        outputter.setOptions(options.output);
        result.statistics = OptimumStatistics();
    }

    void checkInitialGuess()
    {
        if(not lower.empty()) Assert(arma::all(x > lower), "Initial guess is not an interior-point.");
        if(not upper.empty()) Assert(arma::all(x < upper), "Initial guess is not an interior-point.");
        Assert(zl.min() >= 0.0, "Initial guess is not an interior-point.");
        Assert(zu.min() >= 0.0, "Initial guess is not an interior-point.");
    }

    void initialize(double mu_)
    {
        mu = mu_;

        theta0 = arma::norm(h.func, "inf");
        theta_min = 1.0e-4 * std::max(1.0, theta0);
        theta_max = 1.0e+4 * std::max(1.0, theta0);

        clear(filter);
        extend(filter, {theta_max, infinity()});
    }

    void computeNewtonDirection()
    {
        phi.func = f.func - mu * arma::sum(arma::log(x));
        phi.grad = f.grad - mu/x;

        saddle_point_problem.A = h.grad;
        saddle_point_problem.H = f.hessian + arma::diagmat(zl/x);
        saddle_point_problem.f = -(phi.grad - h.grad.t()*y);
        saddle_point_problem.g = -h.func;

        solveFullspaceDense(saddle_point_problem, saddle_point_result, saddle_point_options);

        dx = saddle_point_result.solution.x;
        dy = saddle_point_result.solution.y;
        dz = -(zl % dx + x % zl - mu)/x;
    }

    bool successfulBacktrackingLineSearch()
    {
        // Auxiliary references to ipopt parameters
        const auto& delta                    = options.ipopt.delta;
        const auto& eta_phi                  = options.ipopt.eta_phi;
        const auto& gamma_alpha              = options.ipopt.gamma_alpha;
        const auto& gamma_phi                = options.ipopt.gamma_phi;
        const auto& gamma_theta              = options.ipopt.gamma_theta;
        const auto& s_phi                    = options.ipopt.s_phi;
        const auto& s_theta                  = options.ipopt.s_theta;
        const auto& second_order_corrections = options.ipopt.soc;
        const auto& tau_min                  = options.ipopt.tau_min;

        // Calculate some auxiliary variables for calculating alpha_min
        theta       = arma::norm(h.func, "inf");
        grad_phi_dx = arma::dot(phi.grad, dx);
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
        const double tau = tau_min;
        alphax = fractionToTheBoundary(x, dx, tau);
        alphaz = fractionToTheBoundary(zl, dz, tau);
        alpha  = alphax;

        // Start the backtracking line-search algorithm
        for(unsigned linesearch_iter = 0; ; ++linesearch_iter)
        {
            // Calculate the trial iterate
            x_trial = x + alpha*dx;

            // Update the objective and constraint states with the trial iterate
            f_trial = problem.objective()(x_trial);
            h_trial = problem.constraint()(x_trial);

            // Update the barrier objective function with the trial iterate
            const double phi_trial = f_trial.func - mu * arma::sum(arma::log(x_trial));

            // Compute theta at the trial iterate
            const double theta_trial = arma::norm(h_trial.func, "inf");

            // Check if the new iterate is acceptable in the filter
            const bool filter_acceptable = acceptable(filter, {theta_trial, phi_trial});

            // Calculate some auxiliary variables for checking sufficiency decrease in the measures \theta and \phi
            const double beta_theta = (1 - gamma_theta) * theta;
            const double beta_phi   = phi.func - gamma_phi * theta;

            // Check if the sufficiency decrease condition is satisfied - the condition (18) of the paper
            sufficient_decrease = lessThan(theta_trial, beta_theta, theta_trial) or lessThan(phi_trial, beta_phi, phi_trial);

            // Check if the switching condition is satisfied - the condition (19) of the paper
            switching_condition = grad_phi_dx < 0 and greaterThan(alpha*pow_phi, delta*pow_theta, alpha*pow_phi);

            // Check if the Armijo condition is satisfied - the condition (20) of the paper
            armijo_condition = lessThan(phi_trial - phi.func, eta_phi*alpha*grad_phi_dx, phi.func);

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
                    extend(filter, {beta_theta, beta_phi});
                return true;
            }

            // Check if the second order corrections can be applied and if it succeeds
            if(linesearch_iter == 0 and theta_trial > theta and second_order_corrections)
                if(successfulSecondOrderCorrection())
                    return true;

            // Decrease the length of the size step alpha
            alpha *= 0.5;

            // Check if the step size is smaller than the minimum
            if(alpha < alpha_min)
                return false;
        }
    }

    bool successfulSecondOrderCorrection()
    {
        // Auxiliary references to ipopt parameters
        const auto& max_iters_soc  = options.ipopt.max_iters_soc;
        const auto& tau_min        = options.ipopt.tau_min;
        const auto& kappa_soc      = options.ipopt.kappa_soc;
        const auto& gamma_theta    = options.ipopt.gamma_theta;
        const auto& gamma_phi      = options.ipopt.gamma_phi;
        const auto& eta_phi        = options.ipopt.eta_phi;

        hsoc = alpha*h.func + h_trial.func;

        double theta_soc_old = theta;

        for(unsigned soc_iter = 0; soc_iter < max_iters_soc; ++soc_iter)
        {
            outputter.outputMessage("...applying the second-order correction step");
            saddle_point_problem.g = -hsoc;
            solveFullspaceDense(saddle_point_problem, saddle_point_result, saddle_point_options);
            dx_cor = saddle_point_result.solution.x;
            dy_cor = saddle_point_result.solution.y;

            const double tau = tau_min;
            const double alpha_soc = fractionToTheBoundary(x, dx_cor, tau);

            x_soc = x + alpha_soc * dx_cor;

            f_trial = problem.objective()(x_soc);
            h_trial = problem.constraint()(x_soc);

            // Compute the second-order corrected \theta and \phi measures at the trial iterate
            const double theta_soc = arma::norm(h_trial.func, "inf");
            const double phi_soc = f_trial.func - mu * arma::sum(arma::log(x_soc));

            // Leave the second-order correction algorithm if \theta is not decreasing sufficiently
            if(soc_iter > 0 and theta_soc > kappa_soc * theta_soc_old)
                break;

            // Leave the second-order correction algorithm if x_soc is not acceptable in the filter
            if(not acceptable(filter, {theta_soc, phi_soc}))
                break;

            // Check if the new iterate is acceptable to the filter
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

                outputter.outputMessage("...succeeded!\n");

                // The second order correction was succesfully applied
                return true;
            }

            // Update the constraint evaluation h^soc
            hsoc = alpha_soc*hsoc + h_trial.func;

            // Update the old second-order corrected \theta
            theta_soc_old = theta_soc;
        }

        outputter.outputMessage("...failed!\n");

        // The second order correction could not be succesfully applied
        return false;
    }

    void acceptTrialIterate()
    {
        // Update the next iterate
        x   = x_trial;
        y  +=  alpha * dy;
        zl += alphaz * dz;

        // Update the objective and constraint states
        f = f_trial;
        h = h_trial;
    }

    void updateErrors()
    {
        // Calculate the optimality, feasibility and centrality errors
        errorf = arma::norm(f.grad - h.grad.t()*y - zl, "inf");
        errorh = arma::norm(h.func, "inf");
        errorc = arma::norm((x % zl)/mu - 1, "inf");

        // Calculate the maximum error
        statistics.error = std::max({errorf, errorh, errorc});
    }

    bool converged()
    {
        // Calculte the tolerance for the convergence checking
        const double tol = std::max(options.tolerance, options.ipopt.kappa_epsilon*mu);

        return statistics.error < tol;
    }

    void solveSubproblem(double mu)
    {
        initialize(mu);

        do {
            updateErrors();
            outputState();
            if(converged()) break;
            computeNewtonDirection();
            successfulBacktrackingLineSearch();
            acceptTrialIterate();
            ++statistics.num_iterations;
        } while(statistics.num_iterations < options.max_iterations);

        statistics.converged = statistics.num_iterations < options.max_iterations;
    }

    void minimise()
    {
        f = problem.objective()(x);
        h = problem.constraint()(x);

        outputHeader();

        for(double mu : options.ipopt.mu)
        {
            solveSubproblem(mu);

            if(not statistics.converged) break;
        }
    }

    void outputHeader()
    {
        outputter.addEntry("iter");
        outputter.addEntries("x", n);
        outputter.addEntries("y", m);
        outputter.addEntries("z", n);
        outputter.addEntry("mu");
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
    }

    void outputState()
    {
        outputter.addValue(statistics.num_iterations);
        outputter.addValues(x);
        outputter.addValues(y);
        outputter.addValues(zl);
        outputter.addValue(mu);
        outputter.addValue(f.func);
        outputter.addValue(arma::norm(h.func, "inf"));
        outputter.addValue(errorf);
        outputter.addValue(errorh);
        outputter.addValue(errorc);
        outputter.addValue(statistics.error);
        outputter.addValue(alpha);
        outputter.addValue(alphax);
        outputter.addValue(alphaz);
        outputter.outputState();
    }
};

} // namespace

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
    IpoptSolver solver(problem, result, options);

    solver.minimise();
}

} // namespace Reaktor
