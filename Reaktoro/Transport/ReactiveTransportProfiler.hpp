// // Reaktoro is a unified framework for modeling chemically reactive systems.
// //
// // This library is free software; you can redistribute it and/or
// // modify it under the terms of the GNU Lesser General Public
// // License as published by the Free Software Foundation; either
// // version 2.1 of the License, or (at your option) any later version.
// //
// // This library is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// // Lesser General Public License for more details.
// //
// // You should have received a copy of the GNU Lesser General Public License
// // along with this library. If not, see <http://www.gnu.org/licenses/>.

// #pragma once

// // C++ includes
// #include <ostream>
// #include <memory>

// namespace Reaktoro {

// // Forward declarations
// class ReactiveTransportResult;
// class ReactiveTransportSolver;

// struct ReactiveTransportTimingAccumulated
// {
//     /// The accumulated times for the operations during reactive transport calculations.
//     ReactiveTransportTiming timing_reactive_transport;

//     /// The accumulated times for the operations during equilibrium calculations.
//     EquilibriumTiming timing_equilibrium;

//     /// The accumulated times for the operations during smart equilibrium calculations.
//     SmartEquilibriumTiming timing_smart_equilibrium;
// };

// struct ReactiveTransportProfilingSummary
// {
//     /// The accumulated times for the operations during reactive transport calculations.
//     ReactiveTransportTiming timing_reactive_transport;

//     /// The accumulated times for the operations during equilibrium calculations.
//     EquilibriumTiming timing_equilibrium;

//     /// The accumulated times for the operations during smart equilibrium calculations.
//     SmartEquilibriumTiming timing_smart_equilibrium;

//     /// The total number of chemical equilibrium calculations.
//     Index num_equilibrium_calculations = 0;

//     /// The total number of full chemical equilibrium calculations.
//     Index num_full_equilibrium_calculations = 0;

//     /// The total number of successful smart chemical equilibrium calculations.
//     Index num_smart_equilibrium_calculations = 0;
// };

// struct SmartEquilibriumSuccessTable
// {
//     std::vector<std::vector<bool>> table;
// };

// /// Provide mechanisms for postprocessing and analysis of accumulated profiling/results of a reactive transport simulation.
// class ReactiveTransportProfiler
// {
// public:
//     /// Construct a default instance of ReactiveTransportProfiler.
//     ReactiveTransportProfiler(const ReactiveTransportSolver& solver);

//     /// Construct a copy of a ReactiveTransportProfiler instance.
//     ReactiveTransportProfiler(const ReactiveTransportProfiler& other);

//     /// Destroy this ReactiveTransportProfiler instance.
//     virtual ~ReactiveTransportProfiler();

//     /// Assign a copy of an ReactiveTransportProfiler instance.
//     auto operator=(ReactiveTransportProfiler other) -> ReactiveTransportProfiler&;

//     /// Update the profiler with a new reactive transport time step profiling data.
//     auto update() -> void;

//     /// Return all collected results of the reactive transport calculations.
//     auto results() const -> const std::deque<ReactiveTransportResult>&

//     /// Return a summary of the profiling of the reactive transport calculations.
//     auto summary() const -> ReactiveTransportProfilingSummary;

// private:
//     struct Impl;

//     std::unique_ptr<Impl> pimpl;
// };

// /// Output a summary of the profiling info collected for a reactive transport simulation.
// auto operator<<(std::ostream& out, const ReactiveTransportProfiler& profiler);

// //     /// Name of the file and folder with a status output
// //     std::string folder;
// //     std::string file;

// //     /// Used to store statistics information about the smart equilibrium algorithm.
// //     struct Statistics
// //     {
// //         /// Total time for search operations
// //         double time_estimate = 0.0;

// //         /// Time for search operations (part of estimation)
// //         double time_search = 0.0;

// //         /// Time for matrix-vector multiplications (part of estimation)
// //         double time_mat_vect_mult = 0.0;

// //         /// Time for acceptance test (part of estimation)
// //         double time_acceptance = 0.0;

// //         /// Total time for learn operations
// //         double time_learn = 0.0;

// //         /// Time for store operations (part of learning)
// //         double time_store = 0.0;

// //         /// Time for search operations (part of learning)
// //         double time_gibbs_min = 0.0;

// //         /// The size of the search tree
// //         Index tree_size = 0;

// //         /// Counter for the smart statuses
// //         int learning_counter = 0;

// //         /// Counter for the smart statuses
// //         int total_counter = 0;

// //     };

// //     /// The vector of the  CPU times for learning and estimating time
// //     Statistics total_stats;

// //     /// The vector of the  CPU times for learning and estimating time
// //     std::vector<Statistics> step_stats;

// //     /// The vector of the  CPU times for learning and estimating time
// //     std::vector<bool> statuses;

// //     /// Vector of reactive transport times and equilibrium times
// //     std::vector<double> rt_times;
// //     std::vector<double> eq_times;

// //     /// Flag whether smart of conventional solver was used
// //     bool smart;
// // };

// } // namespace Reaktoro
