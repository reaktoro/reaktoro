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

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>

namespace Reaktoro {

/// Used to provide timing information of the operations during a smart chemical equilibrium calculation.
struct SmartEquilibriumTiming
{
    /// The time spent for solving the chemical equilibrium problem (in seconds).
    double solve = 0.0;

    /// The time spent for learning a new chemical equilibrium state (in seconds).
    double learning = 0.0;

    /// The time spent for a conventional iterative chemical equilibrium calculation during the learning operation (in seconds).
    double learning_solve = 0.0;

    /// The time spent for computing the chemical properties during the learning operation (in seconds).
    double learning_chemical_properties = 0.0;

    /// The time spent for computing the sensitivity matrix during the learning operation (in seconds).
    double learning_sensitivity_matrix = 0.0;

    /// The time spent for computing the error control matrices during the learning operation (in seconds).
    double learning_error_control_matrices = 0.0;

    /// The time spent for storing the computed chemical state into the tree of knowledge (in seconds).
    double learning_storage = 0.0;

    /// The time spent for the smart chemical equilibrium state prediction (in seconds).
    double prediction = 0.0;

    /// The time spent for the search operation during a smart prediction (in seconds).
    double prediction_search = 0.0;

    /// The time spent during on error control while searching during a smart prediction (in seconds).
    double prediction_error_control = 0.0;

    /// The time spent for the matrix-vector multiplication during a smart prediction (in seconds).
    double prediction_taylor = 0.0;

    /// The time spent for updating the priority related info of the clusters in the database (in seconds).
    double prediction_priority_update = 0.0;

    /// Self addition of another SmartEquilibriumTiming instance to this one.
    auto operator+=(const SmartEquilibriumTiming& other) -> SmartEquilibriumTiming&;
};

/// Used to represent the result of a prediction operation in a smart chemical equilibrium calculation.
/// @see SmartEquilibriumResult
struct SmartEquilibriumResultDuringPrediction
{
    /// The indication whether the smart equilibrium prediction was accepted.
    bool accepted = false;

    /// The name of the species that caused the smart approximation to fail.
    String failed_with_species;

    /// The amount of the species that caused the smart approximation to fail.
    double failed_with_amount;

    /// The amount of the species that caused the smart approximation to fail.
    double failed_with_chemical_potential;

    // Self addition assignment to accumulate results.
    auto operator+=(const SmartEquilibriumResultDuringPrediction& other) -> SmartEquilibriumResultDuringPrediction&;
};

/// Used to represent the result of a learning operation in a smart chemical equilibrium calculation.
/// @see SmartEquilibriumResult
struct SmartEquilibriumResultDuringLearning
{
    /// The result of the conventional iterative chemical equilibrium calculation in the learning operation.
    EquilibriumResult solve;

    /// Self addition assignment to accumulate results.
    auto operator+=(const SmartEquilibriumResultDuringLearning& other) -> SmartEquilibriumResultDuringLearning&;
};

/// Used to describe the result of a smart chemical equilibrium calculation.
struct SmartEquilibriumResult
{
    /// Return true if the calculation succeeded.
    auto succeeded() { return prediction.accepted ? true : learning.solve.succeeded(); };

    /// Return true if the calculation failed.
    auto failed() { return !succeeded(); };

    /// Return true if the calculation was performed using a fast first-order Taylor prediction.
    auto predicted() { return prediction.accepted; };

    /// Return true if the calculation was learned, not predicted, and performed using the conventional algorithm.
    auto learned() { return !prediction.accepted; };

    /// Return the number of iterations in the calculation (zero if prediction was successful).
    auto iterations() { return prediction.accepted ? 0 : learning.solve.iterations(); };

    /// The result of the smart approximation operation.
    SmartEquilibriumResultDuringPrediction prediction;

    /// The result of the learning operation (if there was learning).
    SmartEquilibriumResultDuringLearning learning;

    /// The timing information of the operations during a smart chemical equilibrium calculation.
    SmartEquilibriumTiming timing;

    /// Self addition assignment to accumulate results.
    auto operator+=(const SmartEquilibriumResult& other) -> SmartEquilibriumResult&;
};

} // namespace Reaktoro
