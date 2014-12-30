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
#include <Reaktor/Optimization/Filter.hpp>
#include <Reaktor/Optimization/KktSolver.hpp>
#include <Reaktor/Optimization/OptimumProblem.hpp>
#include <Reaktor/Optimization/OptimumOptions.hpp>
#include <Reaktor/Optimization/OptimumResult.hpp>
#include <Reaktor/Optimization/Utils.hpp>

namespace Reaktor {

struct OptimumSolverIpopt::Impl
{
    /// The pointer to the optimisation problem instance
    const OptimumProblem* problem = nullptr;

    /// The pointer to the optimisation result instance
    OptimumResult* result = nullptr;

    /// The pointer to the optimisation options instance
    const OptimumOptions* options = nullptr;

    /// The optimisation filter
    Filter filter;

    /// The components for the the solution of the KKT equations
    KktProblem kkt_problem;
    KktResult kkt_result;
    KktSolver kkt_solver;

    /// The outputter instance
    Outputter outputter;

    /// The auxiliary objective and barrier function results
    ObjectiveResult f, f_trial, phi;

    /// The auxiliary constraint function results
    ConstraintResult h, h_trial;

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

    /// The number of primal variables `n` and equality constraints `m`
    unsigned n, m;

    /// The current value of the barrier parameter
    double mu;

    /// Several other auxiliary algorithmic variables computed along the way
    Vector dl, du, d;
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

    Impl();

    auto initialise(const OptimumProblem& problem, OptimumResult& result, const OptimumOptions& options) -> void;

    auto objective(const Vector& x) -> ObjectiveResult;

    auto constraint(const Vector& x) -> ConstraintResult;

    auto checkInitialGuess() -> void;

    auto restart(double mu_) -> void;

    auto computeNewtonDirection() -> void;

    auto successfulBacktrackingLineSearch() -> bool;

    auto successfulSecondOrderCorrection() -> bool;

    auto acceptTrialIterate() -> void;

    auto updateErrors() -> void;

    auto converged() -> bool;

    auto solveSubproblem(double mu) -> void;

    auto solve(const OptimumProblem& problem, OptimumResult& result, const OptimumOptions& options) -> void;

    auto outputHeader() -> void;

    auto outputState() -> void;
};

OptimumSolverIpopt::Impl::Impl()
{}

auto OptimumSolverIpopt::Impl::initialise(const OptimumProblem& problem, OptimumResult& result, const OptimumOptions& options) -> void
{
    // Set the pointers as the first thing to avoid any segmentation fault
    this->problem = &problem;
    this->result = &result;
    this->options = &options;

    // Set the number of variables `n` and equality constraints `m`
    n = problem.numVariables();
    m = problem.numConstraints();

    // Set the options for the outputter instance
    outputter.setOptions(options.output);

    // Reset the statistics data
    result.statistics = {};

    // Auxiliary references
    auto& x  = result.solution.x;
    auto& zl = result.solution.zl;
    auto& lower = problem.lowerBounds();

    // todo improve this robustness correction against bad initial guess
    x = max(x, lower + options.ipopt.mux*options.ipopt.mu[0]*ones(n));
    zl = options.ipopt.mu[0]/(x - lower);

    // Update the objective and constraint results
    f = objective(x);
    h = constraint(x);
}

auto OptimumSolverIpopt::Impl::objective(const Vector& x) -> ObjectiveResult
{
    Time begin = time();
    ObjectiveResult res = problem->objective()(x);
    result->statistics.time_objective_evals += elapsed(begin);
    return res;
}

auto OptimumSolverIpopt::Impl::constraint(const Vector& x) -> ConstraintResult
{
    Time begin = time();
    ConstraintResult res = problem->constraint()(x);
    result->statistics.time_constraint_evals += elapsed(begin);
    return res;
}

auto OptimumSolverIpopt::Impl::checkInitialGuess() -> void
{
    // Auxiliary references
    const Vector& x  = result->solution.x;
    const Vector& zl = result->solution.zl;
    const Vector& zu = result->solution.zu;
    const Vector& lower = problem->lowerBounds();
    const Vector& upper = problem->upperBounds();

    // TODO improve this annoying way of handling bad initial guesses
    if(lower.size()) if(min(x - lower) > 0.0) error("Cannot continue with the ipopt algorithm", "Initial guess is not an interior-point.");
    if(upper.size()) if(min(upper - x) > 0.0) error("Cannot continue with the ipopt algorithm", "Initial guess is not an interior-point.");
    if(min(zl) >= 0.0) error("Cannot continue with the ipopt algorithm", "Initial guess is not an interior-point.");
    if(min(zu) >= 0.0) error("Cannot continue with the ipopt algorithm", "Initial guess is not an interior-point.");
}

auto OptimumSolverIpopt::Impl::restart(double mu) -> void
{
    this->mu = mu;

    theta0 = norminf(h.func);
    theta_min = 1.0e-4 * std::max(1.0, theta0);
    theta_max = 1.0e+4 * std::max(1.0, theta0);

    filter.clear();
    filter.extend({theta_max, infinity()});
}

auto OptimumSolverIpopt::Impl::computeNewtonDirection() -> void
{
    // Auxiliary references
    const auto& x  = result->solution.x;
    const auto& y  = result->solution.y;
    const auto& zl = result->solution.zl;

    phi.func = f.func - mu * sum(log(x));
    phi.grad = f.grad - mu * inv(x);

    kkt_problem.A = h.grad;
    kkt_problem.H = f.hessian;
    kkt_problem.H.diagonal() += zl/x;
    kkt_problem.f = -(phi.grad - h.grad.transpose()*y);
    kkt_problem.g = -h.func;

    KktOptions kkt_options = options->ipopt.kkt;

    if(options->ipopt.scaling)
        kkt_options.xscaling = sqrt(x);

    kkt_solver.solve(kkt_problem, kkt_result, kkt_options);
    result->statistics.time_linear_system_solutions += kkt_result.statistics.time;

    dx = kkt_result.solution.x;
    dy = kkt_result.solution.y;
    dz = -(zl % dx + x % zl - mu)/x;
}

auto OptimumSolverIpopt::Impl::successfulBacktrackingLineSearch() -> bool
{
    // Auxiliary references
    const auto& x  = result->solution.x;
    const auto& zl = result->solution.zl;

    // Auxiliary references to ipopt parameters
    const auto& delta       = options->ipopt.delta;
    const auto& eta_phi     = options->ipopt.eta_phi;
    const auto& gamma_alpha = options->ipopt.gamma_alpha;
    const auto& gamma_phi   = options->ipopt.gamma_phi;
    const auto& gamma_theta = options->ipopt.gamma_theta;
    const auto& s_phi       = options->ipopt.s_phi;
    const auto& s_theta     = options->ipopt.s_theta;
    const auto& soc         = options->ipopt.soc;
    const auto& tau_min     = options->ipopt.tau_min;

    // Calculate some auxiliary variables for calculating alpha_min
    theta       = norminf(h.func);
    grad_phi_dx = dot(phi.grad, dx);
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
        f_trial = objective(x_trial);
        h_trial = constraint(x_trial);

        // Update the barrier objective function with the trial iterate
        const double phi_trial = f_trial.func - mu * sum(log(x_trial));

        // Compute theta at the trial iterate
        const double theta_trial = norminf(h_trial.func);

        // Check if the new iterate is acceptable in the filter
        const bool filter_acceptable = filter.acceptable({theta_trial, phi_trial});

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
                filter.extend({beta_theta, beta_phi});
            return true;
        }

        // Check if the second order corrections can be applied and if it succeeds
        if(linesearch_iter == 0 and theta_trial > theta and soc)
            if(successfulSecondOrderCorrection())
                return true;

        // Decrease the length of the size step alpha
        alpha *= 0.5;

        // Check if the step size is smaller than the minimum
        if(alpha < alpha_min)
            return false;
    }
}

auto OptimumSolverIpopt::Impl::successfulSecondOrderCorrection() -> bool
{
    // Auxiliary references
    const auto& x = result->solution.x;

    // Auxiliary references to ipopt parameters
    const auto& max_iters_soc  = options->ipopt.max_iters_soc;
    const auto& tau_min        = options->ipopt.tau_min;
    const auto& kappa_soc      = options->ipopt.kappa_soc;
    const auto& gamma_theta    = options->ipopt.gamma_theta;
    const auto& gamma_phi      = options->ipopt.gamma_phi;
    const auto& eta_phi        = options->ipopt.eta_phi;

    h_soc = alpha*h.func + h_trial.func;

    double theta_soc_old = theta;

    for(unsigned soc_iter = 0; soc_iter < max_iters_soc; ++soc_iter)
    {
        outputter.outputMessage("...applying the second-order correction step");
        kkt_problem.g = -h_soc;
        kkt_solver.solve(kkt_problem, kkt_result, options->ipopt.kkt);
        result->statistics.time_linear_system_solutions += kkt_result.statistics.time;
        dx_cor = kkt_result.solution.x;
        dy_cor = kkt_result.solution.y;

        const double tau = tau_min;
        const double alpha_soc = fractionToTheBoundary(x, dx_cor, tau);

        x_soc = x + alpha_soc * dx_cor;

        f_trial = objective(x_soc);
        h_trial = constraint(x_soc);

        // Compute the second-order corrected \theta and \phi measures at the trial iterate
        const double theta_soc = norminf(h_trial.func);
        const double phi_soc = f_trial.func - mu * sum(log(x_soc));

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
        h_soc = alpha_soc * h_soc + h_trial.func;

        // Update the old second-order corrected \theta
        theta_soc_old = theta_soc;
    }

    outputter.outputMessage("...failed!\n");

    // The second order correction could not be succesfully applied
    return false;
}

auto OptimumSolverIpopt::Impl::acceptTrialIterate() -> void
{
    // Update the next iterate
    result->solution.x   = x_trial;
    result->solution.y  += alpha * dy;
    result->solution.zl += alphaz * dz;

    // Update the objective and constraint states
    f = f_trial;
    h = h_trial;
}

auto OptimumSolverIpopt::Impl::updateErrors() -> void
{
    // Auxiliary references
    const auto& x  = result->solution.x;
    const auto& y  = result->solution.y;
    const auto& zl = result->solution.zl;

    // Calculate the optimality, feasibility and centrality errors
    errorf = norminf(f.grad - h.grad.transpose()*y - zl);
    errorh = norminf(h.func);
    errorc = norminf((x % zl)/mu - 1);

    // Calculate the maximum error
    result->statistics.error = std::max({errorf, errorh, errorc});
}

auto OptimumSolverIpopt::Impl::converged() -> bool
{
    // Auxiliary references of the ipopt parameters
    const auto& tolerance = options->tolerance;
    const auto& kappa_epsilon = options->ipopt.kappa_epsilon;

    // Calculte the tolerance for the convergence checking
    const double tol = std::max(tolerance, kappa_epsilon*mu);

    return result->statistics.error < tol;
}

auto OptimumSolverIpopt::Impl::solveSubproblem(double mu) -> void
{
    restart(mu);

    do {
        updateErrors();
        outputState();
        if(converged()) break;
        computeNewtonDirection();
        successfulBacktrackingLineSearch();
        acceptTrialIterate();
        ++result->statistics.num_iterations;
    } while(result->statistics.num_iterations < options->max_iterations);

    // Check if the solution of the subproblem with fixed mu converged
    result->statistics.converged = result->statistics.num_iterations < options->max_iterations;
}

auto OptimumSolverIpopt::Impl::solve(const OptimumProblem& problem, OptimumResult& result, const OptimumOptions& options) -> void
{
    Time begin = time();

    initialise(problem, result, options);

    outputHeader();

    for(double mu : options.ipopt.mu)
    {
        solveSubproblem(mu);

        if(not result.statistics.converged) break;
    }

    result.statistics.time = elapsed(begin);
}

auto OptimumSolverIpopt::Impl::outputHeader() -> void
{
    if(options->output.active)
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
}

auto OptimumSolverIpopt::Impl::outputState() -> void
{
    if(options->output.active)
    {
        // Auxiliary references
        const auto& x  = result->solution.x;
        const auto& y  = result->solution.y;
        const auto& zl = result->solution.zl;

        outputter.addValue(result->statistics.num_iterations);
        outputter.addValues(x);
        outputter.addValues(y);
        outputter.addValues(zl);
        outputter.addValue(mu);
        outputter.addValue(f.func);
        outputter.addValue(norminf(h.func));
        outputter.addValue(errorf);
        outputter.addValue(errorh);
        outputter.addValue(errorc);
        outputter.addValue(result->statistics.error);
        outputter.addValue(alpha);
        outputter.addValue(alphax);
        outputter.addValue(alphaz);
        outputter.outputState();
    }
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
    OptimumOptions default_options;
    pimpl->solve(problem, result, default_options);
}

auto OptimumSolverIpopt::solve(const OptimumProblem& problem, OptimumResult& result, const OptimumOptions& options) -> void
{
    pimpl->solve(problem, result, options);
}


} // namespace Reaktor
