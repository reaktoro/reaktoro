// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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
#include <chrono>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Optimization/OptimumResult.hpp>
#include <Reaktoro/Common/Index.hpp>
namespace Reaktoro {

/// A wrapper class of chrono library to CPU time tracking
struct Timer{

    // Declare the alias for the declare type
    using time_point = std::chrono::time_point<std::chrono::high_resolution_clock>;
    using clock = std::chrono::high_resolution_clock;
    using duration = std::chrono::duration<double>;

    time_point start;
    duration elapsed_time;

    auto startTimer() -> void;
    auto stopTimer() -> double;
};

/// A type used to describe the result of a smart equilibrium calculation.
struct SmartEquilibriumResult
{
    /// The boolean flag that indicates if smart equilibrium calculation was used.
    bool succeeded = false;

    /// Counter for the smart statuses
    std::vector<Index> learning_states_indx;

    /// Used to store statistics information about the smart equilibrium algorithm.
    struct LearnStatistics
    {
        /// Time for learning new state
        double time_learn = 0.0;

        /// Time for reconstructing reference state
        double time_gibbs_min = 0.0;

        /// Time for storing reference state
        double time_store = 0.0;

        /// Apply an addition assignment to this instance
        auto operator+=(const LearnStatistics& other) -> LearnStatistics&;
    };
    struct EstimateStatistics
    {
        /// Time for estimating new state
        double time_estimate = 0.0;

        /// Time for search operations
        double time_search = 0.0;

        /// Time for matrix-vector multiplications
        double time_mat_vect_mult = 0.0;

        /// Time for acceptance test
        double time_acceptance = 0.0;

        /// Apply an addition assignment to this instance
        auto operator+=(const EstimateStatistics& other) -> EstimateStatistics&;

    };

    /// The statistics information of the smart equilibrium algorithm.
    LearnStatistics learn_stats;
    EstimateStatistics estimate_stats;

    /// Apply an addition assignment to this instance
    auto operator+=(const SmartEquilibriumResult& other) -> SmartEquilibriumResult&;
    auto addLearningIndex(const Index & index) -> void;

    /// Counter for the smart statuses
    int smart_counter = 0;

    // The height of the tree storing the reference states
    int tree_size = 0;

};

/// A type used to describe the result of an equilibrium calculation
/// @see ChemicalState
struct EquilibriumResult {
    /// The result of the optimisation calculation
    OptimumResult optimum;

    /// The boolean flag that indicates if smart equilibrium calculation was used.
    SmartEquilibriumResult smart;

    struct Statistics{
        /// Time for learning new state
        double time_learn = 0.0;
        auto operator+=(const Statistics& other) -> Statistics& {
            time_learn     += other.time_learn;
            return *this;
        }
    };
    Statistics stats;

    /// Apply an addition assignment to this instance
    auto operator+=(const EquilibriumResult& other) -> EquilibriumResult&;
};

} // namespace Reaktoro
