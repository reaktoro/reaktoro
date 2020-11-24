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

// Reaktoro includes
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>

namespace Reaktoro {

/// Provide timing information of the operations during a smart equilibrium calculation.
struct SmartEquilibriumTiming
{
    /// The time spent in seconds for solving the chemical equilibrium problem.
    double solve = 0.0;

    /// The time spent in seconds for learning a new chemical equilibrium calculation.
    double learn = 0.0;

    /// The time spent in seconds for a conventional Gibbs energy minimization calculation during learning operation.
    double learn_gibbs_energy_minimization = 0.0;

    /// The time spent in seconds for computing the chemical properties during learning operation.
    double learn_chemical_properties = 0.0;

    /// The time spent in seconds for computing the sensitivity matrix during learning operation.
    double learn_sensitivity_matrix = 0.0;

    /// The time spent in seconds for computing the error control matrices during learning operation.
    double learn_error_control_matrices = 0.0;

    /// The time spent in seconds for storing the computed chemical state into the tree of knowledge.
    double learn_storage = 0.0;

    /// The time spent in seconds for the smart chemical equilibrium state estimation.
    double estimate = 0.0;

    /// The time spent in seconds for the search operation during a smart estimation.
    double estimate_search = 0.0;

    /// The time spent in seconds during on error control while searching during a smart estimation.
    double estimate_error_control = 0.0;

    /// The time spent in seconds for the matrix-vector multiplication during a smart estimation.
    double estimate_taylor = 0.0;

    /// The time spent in seconds for updating the priority related info of the clusters in the database.
    double estimate_database_priority_update = 0.0;

    /// Self addition of another SmartEquilibriumTiming instance to this one.
    auto operator+=(const SmartEquilibriumTiming& other) -> SmartEquilibriumTiming&;
};

/// A type used to define the result status of a smart estimation operation in a smart equilibrium calculation.
/// @see SmartEquilibriumResult
struct SmartEquilibriumResultDuringEstimate
{
    /// The indication whether the smart equilibrium estimate was accepted.
    bool accepted = false;

    /// The name of the species that caused the smart approximation to fail.
    std::string failed_with_species;

    /// The amount of the species that caused the smart approximation to fail.
    double failed_with_amount;

    /// The amount of the species that caused the smart approximation to fail.
    double failed_with_chemical_potential;

    // Self addition assignment to accumulate results.
    auto operator+=(const SmartEquilibriumResultDuringEstimate& other) -> SmartEquilibriumResultDuringEstimate&;
};

/// A type used to define the result status of a learning operation in a smart equilibrium calculation.
/// @see SmartEquilibriumResult
struct SmartEquilibriumResultDuringLearning
{
    /// The result of the full Gibbs energy minimization calculation.
    EquilibriumResult gibbs_energy_minimization;

    /// Self addition assignment to accumulate results.
    auto operator+=(const SmartEquilibriumResultDuringLearning& other) -> SmartEquilibriumResultDuringLearning&;
};

/// A type used to describe the result of a smart equilibrium calculation.
struct SmartEquilibriumResult
{
    /// The result of the smart approximation operation.
    SmartEquilibriumResultDuringEstimate estimate;

    /// The result of the learning operation (if there was learning).
    SmartEquilibriumResultDuringLearning learning;

    /// The timing information of the operations during a smart equilibrium calculation.
    SmartEquilibriumTiming timing;

    /// Self addition assignment to accumulate results.
    auto operator+=(const SmartEquilibriumResult& other) -> SmartEquilibriumResult&;
};

} // namespace Reaktoro
