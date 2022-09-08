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
    learn += other.learn;
    learn_gibbs_energy_minimization += other.learn_gibbs_energy_minimization;
    learn_chemical_properties += other.learn_chemical_properties;
    learn_sensitivity_matrix += other.learn_sensitivity_matrix;
    learn_error_control_matrices += other.learn_error_control_matrices;
    learn_storage += other.learn_storage;
    estimate += other.estimate;
    estimate_search += other.estimate_search;
    estimate_error_control += other.estimate_error_control;
    estimate_taylor += other.estimate_taylor;
    estimate_database_priority_update += other.estimate_database_priority_update;

    return *this;
}

auto SmartEquilibriumResultDuringEstimate::operator+=(const SmartEquilibriumResultDuringEstimate& other) -> SmartEquilibriumResultDuringEstimate&
{
    accepted = other.accepted;
    failed_with_species = other.failed_with_species;
    failed_with_amount = other.failed_with_amount;
    failed_with_chemical_potential = other.failed_with_chemical_potential;

    return *this;
}

auto SmartEquilibriumResultDuringLearning::operator+=(const SmartEquilibriumResultDuringLearning& other) -> SmartEquilibriumResultDuringLearning&
{
    gibbs_energy_minimization +=other.gibbs_energy_minimization;

    return *this;
}

auto SmartEquilibriumResult::operator+=(const SmartEquilibriumResult& other) -> SmartEquilibriumResult&
{
    estimate += other.estimate;
    learning += other.learning;
    timing   += other.timing;

    return *this;
}
} // namespace Reaktoro
