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
#include <Reaktoro/Kinetics/KineticResult.hpp>
#include <Reaktoro/Kinetics/SmartKineticResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>

namespace Reaktoro {

/// Provide a performance analysis of the operations in a kinetic path simulation.
struct KineticPathAnalysis
{
    /// Provide a summary of the performance analysis of all kinetic calculations.
    struct KineticAnalysis
    {
        /// The accumulated timing for the operations during kinetic calculations.
        KineticTiming timing;
    };

    /// Provide a summary of the performance analysis of all smart kinetic calculations.
    struct SmartKineticAnalysis
    {
        /// The accumulated timing for the operations during smart kinetic calculations.
        SmartKineticTiming timing;

        /// The total number of chemical kinetic calculations.
        Index num_kinetic_calculations = 0;

        /// The total number of accepted smart chemical kinetic estimates.
        Index num_smart_kinetic_accepted_estimates = 0;

        /// The total number of required smart chemical kinetic trainings.
        Index num_smart_kinetic_required_learnings = 0;

        /// The success rate at which smart kinetic estimations were accepted.
        double smart_kinetic_estimate_acceptance_rate = 0.0;

        /// The indices of the steps where learning was required.
        std::vector<Index> steps_where_learning_was_required;
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

        /// The indices of the steps where learning was required.
        std::vector<Index> steps_where_learning_was_required;
    };

    /// Provide computing costs (in seconds) for the operations in a kinetic calculation for each time step.
    struct ComputingCostsPerTimeStep
    {
        /// The time spent (in s) in each time step for chemical kinetics calculations.
        std::vector<double> kinetics;

        /// The time spent (in s) in each time step for equilibrium in chemical kinetic calculations.
        std::vector<double> kinetics_equilibration;

        /// The time spent (in s) in each time step for chemical properties evaluation in chemical kinetic calculations.
        std::vector<double> kinetics_properties;

        /// The time spent (in s) in each time step for chemical kinetic calculations without chemical properties evaluation.
        std::vector<double> kinetics_with_ideal_properties;

        /// The time spent (in s) in each time step for smart chemical kinetics calculations.
        std::vector<double> smart_kinetics;

        /// The time spent (in s) in each time step for smart chemical kinetics calculations without the computing costs of search operations.
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

        // The time spent (in s) in each time step for Taylor extrapolation during smart equilibrium estimation operations.
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

    KineticAnalysis kinetics;

    SmartKineticAnalysis smart_kinetics;

    EquilibriumAnalysis equilibrium;

    SmartEquilibriumAnalysis smart_equilibrium;

    ComputingCostsPerTimeStep computing_costs_per_time_step;
};

} // namespace Reaktoro
