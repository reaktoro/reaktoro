// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Optima includes
#include <Optima/Options.hpp>

namespace Reaktoro {

/// The options for the description of the Hessian of the Gibbs energy function
enum class GibbsHessian
{
    /// The Hessian of the Gibbs energy function is fully exact.
    Exact,

    /// The Hessian of the Gibbs energy function is partially exact, partially approximated.
    PartiallyExact,

    /// The Hessian of the Gibbs energy function is approximated using ideal thermodynamic models.
    Approx,

    /// The Hessian of the Gibbs energy function is a diagonal matrix approximation using ideal thermodynamic models.
    ApproxDiagonal,
};

/// The options for the smart equilibrium calculations.
struct SmartEquilibriumOptions
{
    /// The relative tolerance for estimated species mole amounts.
    double reltol = 1.0;

    /// The absolute tolerance for estimated species mole amounts.
    double abstol = 1e-14;
};

/// The options for the equilibrium calculations
struct EquilibriumOptions
{
    /// The options for the optimisation solver.
    Optima::Options optima;

    /// The default lower bound for the amounts of the species.
    double epsilon = 1e-16;

    /// The value multiplied by `epsilon` to compute the logarithm barrier penalty parameter @eq{\tau}.
    double logarithm_barrier_factor = 1.0;

    /// The boolean flag that indicates if warm-start strategy should be used
    /// when possible. Setting this flag to true will cause equilibrium
    /// calculations to use the currect chemical state as an initial guess to
    /// the equilibrium calculation. If the current chemical state is detected
    /// to be uninitiallized (e.g., all species with zero molar amounts), then
    /// cold-start is inevitably used. In this case, a first estimate of a
    /// initial guess will be done using a simplex algorithm, which in most
    /// cases generates a chemical state that works well as initial guess for
    /// all equilibrium algorithms.
    bool warmstart = true;

    /// The calculation mode of the Hessian of the Gibbs energy function
    GibbsHessian hessian = GibbsHessian::PartiallyExact;

    /// The options for the smart equilibrium calculation.
    SmartEquilibriumOptions smart;
};

} // namespace Reaktoro
