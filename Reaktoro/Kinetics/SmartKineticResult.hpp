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

// C++ includes
#include <string>

#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>

namespace Reaktoro {

/// Timing information of the operations during an kinetic calculation.
struct SmartKineticTiming {

    /// The time spent for solving the chemical kinetic problem.
    double solve = 0.0;

    /// The time spent for initializing the chemical kinetic problem.
    double initialize = 0.0;

    /// The time spent for learning the chemical kinetic path.
    double learn = 0.0;

    /// The time spent for integrating the chemical kinetic problem during a smart learning.
    double learn_integration = 0.0;

    /// The time spent for computing the chemical properties of the system during a smart learning.
    double learn_chemical_properties = 0.0;

    /// The time spent for computing the reaction rates of the system during a smart learning.
    double learn_reaction_rates = 0.0;

    /// The time spent for computing the sensitivities of the system during a smart learning.
    double learn_sensitivity = 0.0;

    /// The time spent for computing necessary error matrices a smart learning.
    double learn_error_control_matrices = 0.0;

    /// The time spent for equilibration of the system during a smart learning.
    double learn_equilibration = 0.0;

    /// The time spent for storing the learned smart state.
    double learn_storage = 0.0;

    /// The time spent for estimating the chemical kinetic path.
    double estimate = 0.0;

    /// The time spent for the search operation during a smart estimation.
    double estimate_search = 0.0;

    /// The time spent for the matrix-vector multiplication during a smart estimation.
    double estimate_taylor = 0.0;

    /// The time spent for the acceptance test during a smart estimation.
    double estimate_error_control = 0.0;

    /// The time spent for updating the priority database during a smart estimation.
    double estimate_database_priority_update = 0.0;

    /// The time spent for equilibrating equilibrium species.
    double equilibrate = 0.0;

    /// Self addition assignment to accumulate kinetic timing.
    auto operator+=(const SmartKineticTiming& other) -> SmartKineticTiming&;

};

/// A type used to define the result status of a smart estimation operation in a smart kinetic calculation.
/// @see SmartKineticResult
struct SmartKineticResultDuringEstimate
{
    /// The indication whether the smart equilibrium estimate was accepted.
    bool accepted = false;

    /// The name of the species that caused the smart approximation to fail.
    std::string failed_with_species;

    /// The amount of the species that caused the smart approximation to fail.
    double failed_with_amount;

    /// The amount of the species chemical potential that caused the smart approximation to fail.
    double failed_with_chemical_potential;
};

/// A type used to describe the result of a smart kinetic calculation.
struct SmartKineticResult
{
    /// The timing information of the operations during an equilibrium calculation.
    SmartKineticTiming timing;

    /// The result of the smart approximation operation.
    SmartKineticResultDuringEstimate estimate;

    /// Equilibrium information of the operations during the equilibrium calculation.
    EquilibriumResult equilibrium;

    /// Smart equilibrium information of the operations during smart equilibrium calculation.
    SmartEquilibriumResult smart_equilibrium;

    /// Self addition assignment to accumulate results.
    auto operator+=(const SmartKineticResult& other) -> SmartKineticResult&;
};

} // namespace Reaktoro
