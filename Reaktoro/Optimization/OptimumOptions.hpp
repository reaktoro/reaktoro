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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Outputter.hpp>
#include <Reaktoro/Optimization/KktSolver.hpp>

namespace Reaktoro {

struct OptimumParamsActNewton
{
    /// The threshold for which primal variables lower than it is updated explicitly
    double threshold = 1.0e-14;
};

struct OptimumParamsIpAction
{
    /// The perturbation parameter (or barrier parameter) for the interior-point method
    double mu = 1.0e-16;

    /// The fraction-to-the boundary parameter to relax the line-search backtracking step.
    double tau = 0.9999;

    /// The flag that indicates the optimisation solver to prefer a dense linear system solver even when the Hessian is diagonal
    bool prefer_dense_solver = false;
};

struct OptimumParamsIpNewton
{
    /// The perturbation parameter (or barrier parameter) for the interior-point method
    double mu = 1.0e-8;

    /// The fraction-to-the boundary parameter to relax the line-search backtracking step.
    double tau = 1.0-1e-14;

    /// The factor used to correct the primal initial guess that are too small or on the boundary.
    /// The primal initial guess `x0` is always corrected as `x0' = max(x0, mux*mu)`.
    double mux = 1.0e-5;

    /// The flag that indicates if the KKT problems should be scaled with the matrix sqrt(diag(x))
    bool scaling = true;

    /// The flag that indicates if the direction of the newton step should be used for both primal and dual variables.
    bool uniform_newton_step = false;
};

struct OptimumParamsIpActive
{
    /// The parameter ε for the numerical representation of a zero.
    /// The value of the `i`-th primal variable is considered zero if `x[i] < ε`.
    double epsilon = 1e-20;

    /// The factor τ for the barrier parameter μ defined here as μ = ετ.
    /// The parameter ε is the numerical zero for a primal variable.
    /// @see epsilon
    double tau = 1e-5;
};

struct OptimumParamsIpOpt
{
    /// The interior-point perturbation parameters (default: {1e-8, 1e-16})
    std::vector<double> mu = std::vector<double>({1e-8, 1e-16}); // need explicit initialization list to compile in MSVC 2013

    double delta = 1.0;
    double eta_phi = 1.0e-4;
    double gamma_alpha = 0.05;
    double gamma_phi = 1.0e-5;
    double gamma_theta = 1.0e-5;
    double kappa_epsilon = 10.0;
    double kappa_mu = 0.2;
    double kappa_sigma = 1.0e+10;
    double kappa_soc = 0.99;
    double s_phi = 2.3;
    double s_theta = 1.1;
    double tau_min = 1.0-1e-14;
    double theta_mu = 2.0;
    unsigned max_iters_soc = 4;
    bool soc = true;

    /// The factor used to correct the primal initial guess that are too small or on the boundary.
    double mux = 1.0e-5;

    /// The flag that indicates if the KKT problems should be scaled with the matrix sqrt(diag(x))
    bool scaling = true;
};

struct OptimumParamsKarpov
{
    /// The maximum number of iterations for the line search minimization problem.
    double line_search_max_iterations = 3;

    /// The constant for the Wolfe condition of sufficient decrease in the backtracking line search step
    double line_search_wolfe = 1.0e-4;

    // The fraction-to-the-boundary factor used in the feasible step
    double tau_feasible = 0.99;

    // The fraction-to-the-boundary factor used in the descent step
    double tau_descent = 1.0 - 1.0e-14;

    /// The tolerance for the feasibility problem.
    double feasibility_tolerance = 1.0e-13;

    /// The tolerance for the negative dual variables `z`.
    double negative_dual_tolerance = -1.0e-2;

    /// The value used to remove a primal variable from an active state (i.e., a variable on the bound)
    /// to an inactive state (i.e., a variable in the interior domain).
    double active_to_inactive = 1.0e-6;

    /// The flag that indicates if KktSolver should be used to solve the linear systems
    bool use_kkt_solver = false;
};

struct OptimumParamsRefiner
{
    bool use_lma_setup = true;
};

/// A type that describes the options for the output of a optimisation calculation
struct OptimumOutput : OutputterOptions
{
    /// The prefix for the primal variables `x`.
    std::string xprefix = "x";

    /// The prefix for the dual variables `y`.
    std::string yprefix = "y";

    /// The prefix for the dual variables `z`.
    std::string zprefix = "z";

    /// The names of the primal variables `x`.
    /// Numbers will be used if not properly set (e.g., `x[0]`, `x[1]`)
    std::vector<std::string> xnames;

    /// The names of the dual variables `y`.
    /// Numbers will be used if not properly set (e.g., `y[0]`, `y[1]`)
    std::vector<std::string> ynames;

    /// The names of the dual variables `z`.
    /// Numbers will be used if not properly set (e.g., `z[0]`, `z[1]`)
    std::vector<std::string> znames;
};

/// A type that describes the options of a optimisation calculation
struct OptimumOptions
{
    /// The residual tolerance in the optimisation calculations
    double tolerance = 1.0e-6;

    /// The maximum number of iterations in the optimisation calculations
    unsigned max_iterations = 2000;

    /// The parameters for the ActNewton algorithm
    OptimumParamsActNewton actnewton;

    /// The options for the output of the optimisation calculations
    OptimumOutput output;

    /// The parameters for the IpAction algorithm
    OptimumParamsIpAction ipaction;

    /// The parameters for the IpOpt algorithm
    OptimumParamsIpOpt ipopt;

    /// The parameters for the IpNewton algorithm
    OptimumParamsIpNewton ipnewton;

    /// The parameters for the IpActive algorithm
    OptimumParamsIpActive ipactive;

    /// The parameters for the Karpov algorithm
    OptimumParamsKarpov karpov;

    /// The parameters for the Refiner algorithm
    OptimumParamsRefiner refiner;

    /// The options for the KKT calculations
    KktOptions kkt;
};

} // namespace Reaktoro
