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

#pragma once

// C++ includes
#include <vector>

// Reaktor includes
#include <Reaktor/Optimization/SaddlePointUtils.hpp>

namespace Reaktor {

// Forward declarations
struct OptimumOptions;
class OptimumProblem;
struct OptimumResult;

struct IpoptOptions
{
    std::vector<double> mu = {1e-8, 1e-16};
    double delta           = 1.0;
    double eta_phi         = 1.0e-4;
    double gamma_alpha     = 0.05;
    double gamma_phi       = 1.0e-5;
    double gamma_theta     = 1.0e-5;
    double kappa_epsilon   = 10.0;
    double kappa_mu        = 0.2;
    double kappa_sigma     = 1.0e+10;
    double kappa_soc       = 0.99;
    double s_phi           = 2.3;
    double s_theta         = 1.1;
    double tau_min         = 0.99;
    double theta_mu        = 2.0;
    unsigned max_iters_soc = 4;
    bool soc               = true;

    /// The factor used to correct the primal initial guess that are too small or on the boundary.
    double mux = 1.0e-5;

    /// The options for the saddle point problem calculations
    SaddlePointOptions saddle_point;
};

auto ipopt(const OptimumProblem& problem, OptimumResult& result, const OptimumOptions& options) -> void;

} // namespace Reaktor
