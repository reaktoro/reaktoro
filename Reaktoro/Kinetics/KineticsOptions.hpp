// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

// Reaktoro includes
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>

namespace Reaktoro {

/// The options for chemical kinetics calculation.
struct KineticsOptions : EquilibriumOptions
{
    /// Construct a default KineticsOptions object.
    KineticsOptions()
    {
        // For kinetics, it is crucial that every calculation performs at least
        // one Newton step, even if the 0th iteration already shows sactisfatory
        // residuals that pass the convergence test. Without this, we see
        // simulations in which, for example, a mineral initially absent (with
        // amount 1e-16) that should precipitate ends up not precipitating
        // because the calculation already exists at the 0th iteration, without
        // a Newton step applied.
        EquilibriumOptions::optima.convergence.requires_at_least_one_iteration = true;

        // TODO: Requiring at least one Newton step in EquilibriumSolver should
        // be default. This needs Optima to be adjusted, however, so that it
        // always pick the same set of basic variables that resulted from the
        // last calculation. Currently, this does not happen, and some
        // EquilibriumSolver tests fail because it was expecting just an extra
        // iteration for convergence, but becasue of the change in basic
        // variables, 3 or 4 iterations become needed instead.
    }

    /// Construct a  KineticsOptions object from a EquilibriumOptions one.
    KineticsOptions(EquilibriumOptions const& other)
    : EquilibriumOptions(other)
    {
        // Check reasons above.
        EquilibriumOptions::optima.convergence.requires_at_least_one_iteration = true;
    }

    /// The time step used for preconditioning the chemical state when performing the very first chemical kinetics step.
    double dt0 = 1e-6;
};

} // namespace Reaktoro
