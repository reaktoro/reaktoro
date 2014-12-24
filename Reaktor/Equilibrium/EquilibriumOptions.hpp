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
#include <Reaktor/Optimization/OptimumOptions.hpp>

namespace Reaktor {

/// The numerical optimization algorithm used for the chemical equilibrium calculations
enum OptimizationAlgorithm
{
    IpnewtonAlgorithm, IpoptAlgorithm
};

/// The flag indicating how the second order derivatives of the Gibbs energy function is handled.
/// Currently, only NumericalHessian and DiagonalHessian are supported.
enum HessianApproximation
{
    DiagonalHessian, AnalyticalHessian, IdealHessian
};

/// The options for the equilibrium calculations
struct EquilibriumOptions
{
    /// Construct a default EquilibriumOptions instance
    EquilibriumOptions();

    /// The choice of numerical optimization algorithm for the minimization of the Gibbs energy function
    OptimizationAlgorithm algorithm = IpnewtonAlgorithm;

    /// The choice of hessian approaximation for the Gibbs energy function
    HessianApproximation hessian = DiagonalHessian;

    /// The options for the optimization calculation.
    OptimumOptions optimization;

    /// The boolean flag that indicates if the options in `optimization` should be automatically tunned.
    /// In the *smart mode*, the options in `optimization` will be adjusted accordingly with the
    /// given `algorithm` and 'hessian' options.
    bool smart_mode = true;
};

} // namespace Reaktor
