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
#include <deque>
#include <ostream>
#include <memory>

// Reaktoro includes
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>
#include <Reaktoro/Transport/ReactiveTransportResult.hpp>
#include <Reaktoro/Transport/TransportResult.hpp>

namespace Reaktoro {

// Forward declarations
class ReactiveTransportSolver;

/// Provide mechanisms for analysing accumulated profiling/results of a reactive transport calculation.
class ReactiveTransportProfiler
{
public:
    /// Provide accumulated timing information in each time step for the operations in a reactive transport calculation.
    struct ComputingCostsPerTimeStep
    {
        /// The time (in s) for each time step.
        std::vector<double> t;

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

        /// The time spent (in s) in each time step for Gibbs energy minimization calculations during smart equilibrium learning operations.
        std::vector<double> smart_equilibrium_gibbs_energy_minimization;

        /// The time spent (in s) in each time step for storing a new learned chemical state in smart chemical equilibrium calculations.
        std::vector<double> smart_equilibrium_storage;
    };

    /// Provide a summary of the performance analysis of the smart equilibrium operations in a reactive transport calculation.
    struct SmartEquilibriumProfiling
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

    /// Construct a default instance of ReactiveTransportProfiler.
    ReactiveTransportProfiler(const ReactiveTransportSolver& solver);

    /// Construct a copy of a ReactiveTransportProfiler instance.
    ReactiveTransportProfiler(const ReactiveTransportProfiler& other);

    /// Destroy this ReactiveTransportProfiler instance.
    virtual ~ReactiveTransportProfiler();

    /// Assign a copy of an ReactiveTransportProfiler instance.
    auto operator=(ReactiveTransportProfiler other) -> ReactiveTransportProfiler&;

    /// Update the profiler with a new reactive transport time step profiling data.
    auto update() -> void;

    /// Return the computing costs of all operations during a reactive transport calculation.
    auto computingCostsPerTimeStep() const -> ComputingCostsPerTimeStep;

    /// Return a summary of the performance analysis of the smart equilibrium operations in a reactive transport calculation.
    auto smartEquilibriumProfiling() const -> SmartEquilibriumProfiling;

    /// Output the complete analysis of the performance of reactive transport simulation.
    auto output(std::string filename) -> void;

    /// Return all collected results of the reactive transport calculations.
    auto results() const -> const std::deque<ReactiveTransportResult>&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

/// Output a summary of the profiling info collected for a reactive transport simulation.
auto operator<<(std::ostream& out, const ReactiveTransportProfiler& profiler);

/// Output the collected profiling information of the smart chemical equilibrium calculations in a reactive transport simulation.
auto operator<<(std::ostream& out, const ReactiveTransportProfiler::SmartEquilibriumProfiling& smartprof) -> std::ostream&;

/// Output the table of computing costs per time step of the operations in a reactive transport simulation.
auto operator<<(std::ostream& out, const ReactiveTransportProfiler::ComputingCostsPerTimeStep& costs) -> std::ostream&;




//     /// Name of the file and folder with a status output
//     std::string folder;
//     std::string file;

//     /// Used to store statistics information about the smart equilibrium algorithm.
//     struct Statistics
//     {
//         /// Total time for search operations
//         double time_estimate = 0.0;

//         /// Time for search operations (part of estimation)
//         double time_search = 0.0;

//         /// Time for matrix-vector multiplications (part of estimation)
//         double time_mat_vect_mult = 0.0;

//         /// Time for acceptance test (part of estimation)
//         double time_acceptance = 0.0;

//         /// Total time for learn operations
//         double time_learn = 0.0;

//         /// Time for store operations (part of learning)
//         double time_store = 0.0;

//         /// Time for search operations (part of learning)
//         double time_gibbs_min = 0.0;

//         /// The size of the search tree
//         Index tree_size = 0;

//         /// Counter for the smart statuses
//         int learning_counter = 0;

//         /// Counter for the smart statuses
//         int total_counter = 0;

//     };

//     /// The vector of the  CPU times for learning and estimating time
//     Statistics total_stats;

//     /// The vector of the  CPU times for learning and estimating time
//     std::vector<Statistics> step_stats;

//     /// The vector of the  CPU times for learning and estimating time
//     std::vector<bool> statuses;

//     /// Vector of reactive transport times and equilibrium times
//     std::vector<double> rt_times;
//     std::vector<double> eq_times;

//     /// Flag whether smart of conventional solver was used
//     bool smart;
// };

} // namespace Reaktoro
