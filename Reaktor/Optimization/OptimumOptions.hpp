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
#include <Reaktor/Common/Outputter.hpp>
#include <Reaktor/Optimization/AlgorithmIpnewton.hpp>
#include <Reaktor/Optimization/AlgorithmIpopt.hpp>
#include <Reaktor/Optimization/SaddlePointUtils.hpp>

namespace Reaktor {

/// A type that describes
struct OptimumOptions
{
    /// The residual tolerance in the optimization calculation
    double tolerance = 1.0e-8;

    /// The maximum number of iterations in the optimization calculation
    unsigned max_iterations = 200;

    /// The options for the output of the optimization calculation
    OutputOptions output;

    /// The options for the optimization calculation when using the ipopt algorithm
    IpoptOptions ipopt;

    /// The options for the optimization calculation when using the ipnewton algorithm
    IpnewtonOptions ipnewton;

    /// The options for the saddle point problem calculations
    SaddlePointOptions saddlepoint;
};

} // namespace Reaktor
