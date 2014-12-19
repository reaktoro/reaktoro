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

// Reaktor includes
#include <Reaktor/Optimization/SaddlePointUtils.hpp>

namespace Reaktor {

// Forward declarations
struct OptimumOptions;
class OptimumProblem;
struct OptimumResult;

struct IpnewtonOptions
{
    /// The perturbation parameter (or barrier parameter) for the interior-point method
    double mu = 1.0e-8;

    /// The fraction-to-the boundary parameter to relax the line-search backtracking step
    double tau = 0.99;

    /// The factor used to correct the primal initial guess that are too small or on the boundary.
    /// The primal initial guess `x0` is always corrected as `x0' = max(x0, mux*mu)`.
    double mux = 1.0e-5;

    /// The options for the saddle point problem calculations
    SaddlePointOptions saddle_point;

    /// The flag that indicates if the saddle point problems should be scaled with the matrix sqrt(diag(x))
    bool scaling = true;

    /// The flag that indicates if the direction of the newton step should be used for both primal and dual variables
    bool uniform_newton_step = true;
};

auto ipnewton(const OptimumProblem& problem, OptimumResult& result, const OptimumOptions& options) -> void;

} // namespace Reaktor
