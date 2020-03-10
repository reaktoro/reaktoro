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
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>
#include <Reaktoro/Transport/TransportResult.hpp>

namespace Reaktoro {

/// Provide a performance analysis of the operations in a reactive transport simulation.
struct ReactiveTransportAnalysis
{
    /// Provide a summary of the performance analysis of all transport calculations.
    struct TransportAnalysis
    {
        /// The accumulated timing for the operations during fluid element transport calculations.
        TransportTiming timing;
    };

    /// Provide a summary of the performance analysis of all equilibrium calculations.
    struct EquilibriumAnalysis
    {
        /// The accumulated timing for the operations during equilibrium calculations.
        EquilibriumTiming timing;
    };

    /// Provide a summary of the performance analysis of all smart equilibrium calculations.
    struct SmartEquilibriumAnalysis
    {
        /// The accumulated timing for the operations during smart equilibrium calculations.
        SmartEquilibriumTiming timing;

        /// The total number of chemical equilibrium calculations.
        Index num_equilibrium_calculations = 0;

        /// The total number of accepted smart chemical equilibrium estimates.
        Index num_smart_equilibrium_accepted_estimates = 0;

        /// The total number of required smart chemical equilibrium trainings.
        Index num_smart_equilibrium_required_learnings = 0;

        /// The success rate at which smart equilibrium estimates were accepted.
        double smart_equilibrium_estimate_acceptance_rate = 0.0;

        /// The indices of the cell at each time step where learning was required.
        std::vector<std::vector<Index>> cells_where_learning_was_required_at_step;
    };

    /// Provide computing costs (in seconds) for the operations in a reactive transport calculation for each time step.
    struct ComputingCostsPerTimeStep
    {
        /// The time spent (in s) in each time step for fluid element transport calculations.
        std::vector<double> transport;

        /// The time spent (in s) in each time step for chemical equilibrium calculations.
        std::vector<double> equilibrium;

        /// The time spent (in s) in each time step for smart chemical equilibrium calculations.
        std::vector<double> smart_equilibrium;

        /// The time spent (in s) in each time step for smart chemical equilibrium calculations without the computing costs of nearest neighbor search operations.
        std::vector<double> smart_equilibrium_with_ideal_search;

        /// The time spent (in s) in each time step for smart chemical equilibrium estimation calculations.
        std::vector<double> smart_equilibrium_estimate;

        /// The time spent (in s) in each time step for nearest neighbor search operations during smart chemical equilibrium calculations.
        std::vector<double> smart_equilibrium_nearest_neighbor_search;

        /// The time spent (in s) in each time step for smart chemical equilibrium estimation calculations.
        std::vector<double> smart_equilibrium_database_priority_update;

        /// The time spent (in s) in each time step for Gibbs energy minimization calculations during smart equilibrium learning operations.
        std::vector<double> smart_equilibrium_gibbs_energy_minimization;

        /// The time spent (in s) in each time step for storing a new learned chemical state in smart chemical equilibrium calculations.
        std::vector<double> smart_equilibrium_storage;
    };

    TransportAnalysis transport;

    EquilibriumAnalysis equilibrium;

    SmartEquilibriumAnalysis smart_equilibrium;

    ComputingCostsPerTimeStep computing_costs_per_time_step;
};

} // namespace Reaktoro
