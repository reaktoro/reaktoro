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

#include "SmartEquilibriumResult.hpp"

namespace Reaktoro {

auto SmartEquilibriumTiming::operator+=(const SmartEquilibriumTiming& other) -> SmartEquilibriumTiming&
{
    solve += other.solve;
    learning += other.learning;
    learning_solve += other.learning_solve;
    learning_chemical_properties += other.learning_chemical_properties;
    learning_sensitivity_matrix += other.learning_sensitivity_matrix;
    learning_error_control_matrices += other.learning_error_control_matrices;
    learning_storage += other.learning_storage;
    prediction += other.prediction;
    prediction_search += other.prediction_search;
    prediction_error_control += other.prediction_error_control;
    prediction_taylor += other.prediction_taylor;
    prediction_priority_update += other.prediction_priority_update;

    return *this;
}

auto SmartEquilibriumResultDuringPrediction::operator+=(const SmartEquilibriumResultDuringPrediction& other) -> SmartEquilibriumResultDuringPrediction&
{
    accepted = other.accepted;
    failed_with_species = other.failed_with_species;
    failed_with_amount = other.failed_with_amount;
    failed_with_chemical_potential = other.failed_with_chemical_potential;

    return *this;
}

auto SmartEquilibriumResultDuringLearning::operator+=(const SmartEquilibriumResultDuringLearning& other) -> SmartEquilibriumResultDuringLearning&
{
    solve +=other.solve;

    return *this;
}

auto SmartEquilibriumResult::operator+=(const SmartEquilibriumResult& other) -> SmartEquilibriumResult&
{
    prediction += other.prediction;
    learning += other.learning;
    timing   += other.timing;

    return *this;
}
} // namespace Reaktoro
