// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

    /// Provide computing costs (in seconds) for the operations in a reactive transport calculation for each time step.
    struct ComputingCostsPerTimeStep
    {
        /// The time spent (in s) in each time step for fluid element transport calculations.
        std::vector<double> transport;

        /// The time spent (in s) in each time step for chemical kinetic calculations.
        std::vector<double> kinetics;

        /// The time spent (in s) in each time step for equilibrium in chemical kinetic calculations.
        std::vector<double> kinetics_equilibration;

        /// The time spent (in s) in each time step for chemical properties evaluation in chemical kinetic calculations.
        std::vector<double> kinetics_properties;

        /// The time spent (in s) in each time step for chemical kinetic calculations without chemical properties evaluation.
        std::vector<double> kinetics_with_ideal_properties;

        /// The time spent (in s) in each time step for smart chemical kinetic calculations.
        std::vector<double> smart_kinetics;

        /// The time spent (in s) in each time step for smart chemical kinetics calculations without the computing costs of nearest neighbor search operations.
        std::vector<double> smart_kinetics_with_ideal_search;

        /// The time spent (in s) in each time step for learning in smart chemical kinetic calculations.
        std::vector<double> smart_kinetics_learn;

        /// The time spent (in s) in each time step for chemical properties evaluation as part of learning in smart chemical kinetic calculations.
        std::vector<double> smart_kinetics_chemical_properties;

        /// The time spent (in s) in each time step for chemical properties evaluation as part of learning in smart chemical kinetic calculations.
        std::vector<double> smart_kinetics_equilibration;

        /// The time spent (in s) in each time step for estimating in smart chemical kinetic calculations.
        std::vector<double> smart_kinetics_estimate;

        /// The time spent (in s) in each time step for searching as part of estimating in smart chemical kinetic calculations.
        std::vector<double> smart_kinetics_search;

        /// The time spent (in s) in each time step for acceptance verification as part of estimating in smart chemical kinetic calculations.
        std::vector<double> smart_kinetics_error_control;

        /// The time spent (in s) in each time step for matrix-vector multiplication as part of estimating in smart chemical kinetic calculations.
        std::vector<double> smart_kinetics_taylor;

        /// The time spent (in s) in each time step for chemical equilibrium calculations.
        std::vector<double> equilibrium;

        /// The time spent (in s) in each time step for smart chemical equilibrium calculations.
        std::vector<double> smart_equilibrium;

        /// The time spent (in s) in each time step for smart chemical equilibrium calculations without the computing costs of search operations.
        std::vector<double> smart_equilibrium_with_ideal_search;

        /// The time spent (in s) in each time step for smart chemical equilibrium estimation calculations.
        std::vector<double> smart_equilibrium_estimate;

        /// The time spent (in s) in each time step for search operations during smart equilibrium estimation operations.
        std::vector<double> smart_equilibrium_search;

        /// The time spent (in s) in each time step for error control operations during smart equilibrium estimation operations.
        std::vector<double> smart_equilibrium_error_control;

        /// The time spent (in s) in each time step for Taylor extrapolation during smart equilibrium estimation operations.
        std::vector<double> smart_equilibrium_taylor;

        /// The time spent (in s) in each time step for updating the priority related info during smart equilibrium estimation operations.
        std::vector<double> smart_equilibrium_database_priority_update;

        /// The time spent (in s) in each time step for smart equilibrium learning calculations.
        std::vector<double> smart_equilibrium_learn;

        /// The time spent (in s) in each time step for Gibbs energy minimization calculations during smart equilibrium learning operations.
        std::vector<double> smart_equilibrium_gibbs_energy_minimization;

        /// The time spent (in s) in each time step for computing the chemical properties during smart equilibrium learning operation.
        std::vector<double> smart_equilibrium_chemical_properties;

        /// The time spent (in s) in each time step  for computing the sensitivity matrix during smart equilibrium learning operation.
        std::vector<double> smart_equilibrium_sensitivity_matrix;

        /// The time spent (in s) in each time step for computing the error control matrices during smart equilibrium learning operation.
        std::vector<double> smart_equilibrium_error_control_matrices;

        /// The time spent (in s) in each time step for storing a new learned chemical state during smart equilibrium learning operation.
        std::vector<double> smart_equilibrium_storage;
    };

    TransportAnalysis transport;

    EquilibriumAnalysis equilibrium;

    ComputingCostsPerTimeStep computing_costs_per_time_step;
};

} // namespace Reaktoro
