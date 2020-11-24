// Reaktoro is a unified framework for modeling chemically reactive systems.
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

#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>

namespace Reaktoro {

/// Timing information of the operations during an kinetic calculation.
struct KineticTiming {

    /// The time spent for solving the chemical kinetic problem.
    double solve = 0.0;

    /// The time spent for initializing the chemical kinetic problem.
    double initialize = 0.0;

    /// The time spent for integrating the chemical kinetic problem.
    double integrate = 0.0;

    /// The time spent for computing the chemical properties of the system during the integration step.
    double integrate_chemical_properties = 0.0;

    /// The time spent for computing the reaction rates of the system during the integration step.
    double integrate_reaction_rates = 0.0;

    /// The time spent for computing the sensitivities of the system during the integration step.
    double integrate_sensitivity = 0.0;

    /// The time spent for equilibration of the system triggered during evaluation of right-hand side during the integration step.
    double integrate_equilibration = 0.0;

    /// The time spent for equilibration of the species after integration of kinetic species.
    double equilibrate = 0.0;

    /// Self addition assignment to accumulate kinetic timing.
    auto operator+=(const KineticTiming& other) -> KineticTiming&;

};
/// A type used to describe the result of an kinetic calculation.
struct KineticResult
{
    /// The timing information of the operations during the kinetic calculation.
    KineticTiming timing;

    /// Equilibrium information of the operations during the kinetic calculation.
    EquilibriumResult equilibrium;

    /// Smart equilibrium information of the operations during the kinetic calculation.
    SmartEquilibriumResult smart_equilibrium;

    /// Self addition assignment to accumulate results.
    auto operator+=(const KineticResult& other) -> KineticResult&;
};

} // namespace Reaktoro
