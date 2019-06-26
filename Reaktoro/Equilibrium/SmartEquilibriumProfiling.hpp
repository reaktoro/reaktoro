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

namespace Reaktoro {

/// Provide profiling information of the operations during a smart equilibrium calculation.
struct SmartEquilibriumProfiling
{
    /// The total time spent for learning a new chemical equilibrium calculation.
    double time_learning = 0.0;

    /// The time spent for a conventional Gibbs energy minimization calculation during learning operation.
    double time_learning_gibbs_energy_minimization = 0.0;

    /// The time spent for computing the chemical properties during learning operation.
    double time_learning_chemical_properties = 0.0;

    /// The time spent for computing the sensitivity matrix during learning operation.
    double time_learning_sensitivity_matrix = 0.0;

    /// The time spent for storing the computed chemical state into the tree of knowledge.
    double time_learning_storage = 0.0;

    /// The total time spent for the smart chemical equilibrium state estimation.
    double time_estimate = 0.0;

    /// The time spent for the search operation during a smart estimation.
    double time_estimate_search = 0.0;

    /// The time spent for the matrix-vector multiplication during a smart estimation.
    double time_estimate_mat_vec_mul = 0.0;

    /// The time spent for the acceptance test during a smart estimation.
    double time_estimate_acceptance = 0.0;
};

} // namespace Reaktoro
