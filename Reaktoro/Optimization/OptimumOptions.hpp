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

struct OptimumParamsIpnewton
{
    /// The perturbation parameter (or barrier parameter) for the interior-point method
    double mu = 1.0e-8;

    /// The fraction-to-the boundary parameter to relax the line-search backtracking step.
    double tau = 0.999995;

    /// The factor used to correct the primal initial guess that are too small or on the boundary.
    /// The primal initial guess `x0` is always corrected as `x0' = max(x0, mux*mu)`.
    double mux = 1.0e-5;

    /// The flag that indicates if the KKT problems should be scaled with the matrix sqrt(diag(x))
    bool scaling = true;

    /// The flag that indicates if the direction of the newton step should be used for both primal and dual variables.
    /// Using a uniform Newton step is usually more robust and converges faster for warm start initial guesses.
    bool uniform_newton_step = true;
};

struct OptimumParamsIpopt
{
    std::vector<double> mu = { 1e-8, 1e-16 };

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
    double tau_min = 0.99;
    double theta_mu = 2.0;
    unsigned max_iters_soc = 4;
    bool soc = true;

    /// The factor used to correct the primal initial guess that are too small or on the boundary.
    double mux = 1.0e-5;

    /// The flag that indicates if the KKT problems should be scaled with the matrix sqrt(diag(x))
    bool scaling = true;
};

/// The method used for the optimisation calculationss
enum class OptimumMethod
{
    Ipnewton, Ipopt
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
    double tolerance = 1.0e-8;

    /// The maximum number of iterations in the optimisation calculations
    unsigned max_iterations = 5000;

    /// The algorithm for the optimisation calculations
    OptimumMethod method = OptimumMethod::Ipnewton;

    /// The options for the output of the optimisation calculations
    OptimumOutput output;

    /// The parameters for the ipopt algorithm
    OptimumParamsIpopt ipopt;

    /// The parameters for the ipnewton algorithm
    OptimumParamsIpnewton ipnewton;

    /// The options for the KKT calculations
    KktOptions kkt;
};

} // namespace Reaktoro
