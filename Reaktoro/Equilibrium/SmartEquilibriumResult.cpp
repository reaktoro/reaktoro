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

#include "SmartEquilibriumResult.hpp"

namespace Reaktoro {

auto SmartEquilibriumTiming::operator+=(const SmartEquilibriumTiming& other) -> SmartEquilibriumTiming&
{
    solve += other.solve;
    learn += other.learn;
    learning_gibbs_energy_minimization += other.learning_gibbs_energy_minimization;
    learning_chemical_properties += other.learning_chemical_properties;
    learning_sensitivity_matrix += other.learning_sensitivity_matrix;
    learning_error_control_matrices += other.learning_error_control_matrices;
    learning_storage += other.learning_storage;
    estimate += other.estimate;
    estimate_search += other.estimate_search;
    estimate_error_control += other.estimate_error_control;
    estimate_taylor += other.estimate_taylor;
    estimate_database_priority_update += other.estimate_database_priority_update;
    return *this;
}

} // namespace Reaktoro
