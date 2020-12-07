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

#include "JsonTransport.hpp"

// Reaktoro includes
#include <Reaktoro/Serialization/JsonEquilibrium.hpp>
#include <Reaktoro/Serialization/JsonKinetics.hpp>
#include <Reaktoro/Transport/TransportResult.hpp>
#include <Reaktoro/Transport/ReactiveTransportAnalysis.hpp>

namespace Reaktoro {

void to_json(json& j, const TransportTiming& obj) {
    j["step"] = obj.step;
    j["matrix_equation_assembly"] = obj.matrix_equation_assembly;
    j["matrix_equation_solve"] = obj.matrix_equation_solve;
}

void from_json(const json& j, TransportTiming& obj) {
    j.at("step").get_to(obj.step);
    j.at("matrix_equation_assembly").get_to(obj.matrix_equation_assembly);
    j.at("matrix_equation_solve").get_to(obj.matrix_equation_solve);
}

void to_json(json& j, const TransportResult& obj) {
    j["timing"] = obj.timing;
}

void from_json(const json& j, TransportResult& obj) {
    j.at("timing").get_to(obj.timing);
}

void to_json(json& j, const ReactiveTransportAnalysis& obj) {
    j["transport"] = json::object();
    j["transport"]["timing"] = obj.transport.timing;

    j["equilibrium"] = json::object();
    j["equilibrium"]["timing"] = obj.equilibrium.timing;

    j["smart_equilibrium"] = json::object();
    j["smart_equilibrium"]["timing"] = obj.smart_equilibrium.timing;
    j["num_equilibrium_calculations"] = obj.smart_equilibrium.num_equilibrium_calculations;
    j["num_smart_equilibrium_accepted_estimates"] = obj.smart_equilibrium.num_smart_equilibrium_accepted_estimates;
    j["num_smart_equilibrium_required_learnings"] = obj.smart_equilibrium.num_smart_equilibrium_required_learnings;
    j["smart_equilibrium_estimate_acceptance_rate"] = obj.smart_equilibrium.smart_equilibrium_estimate_acceptance_rate;
    j["smart_equilibrium_cells_where_learning_was_required_at_step"] = obj.smart_equilibrium.cells_where_learning_was_required_at_step;

    j["kinetics"] = json::object();
    j["kinetics"]["timing"] = obj.kinetics.timing;

    j["smart_kinetics"] = json::object();
    j["smart_kinetics"]["timing"] = obj.smart_kinetics.timing;
    j["num_kinetics_calculations"] = obj.smart_kinetics.num_kinetics_calculations;
    j["num_smart_kinetics_accepted_estimates"] = obj.smart_kinetics.num_smart_kinetics_accepted_estimates;
    j["num_smart_kinetics_required_learnings"] = obj.smart_kinetics.num_smart_kinetics_required_learnings;
    j["smart_kinetics_estimate_acceptance_rate"] = obj.smart_kinetics.smart_kinetics_estimate_acceptance_rate;
    j["smart_kinetics_cells_where_learning_was_required_at_step"] = obj.smart_kinetics.cells_where_learning_was_required_at_step;

    j["computing_costs_per_time_step"] = json::object();
    j["computing_costs_per_time_step"]["transport"] = obj.computing_costs_per_time_step.transport;
    j["computing_costs_per_time_step"]["equilibrium"] = obj.computing_costs_per_time_step.equilibrium;
    j["computing_costs_per_time_step"]["kinetics"] = obj.computing_costs_per_time_step.kinetics;

    j["computing_costs_per_time_step"]["smart_kinetics"] = obj.computing_costs_per_time_step.smart_kinetics;
    j["computing_costs_per_time_step"]["smart_kinetics_with_ideal_search"] = obj.computing_costs_per_time_step.smart_kinetics_with_ideal_search;

    j["computing_costs_per_time_step"]["smart_kinetics_estimate"] = obj.computing_costs_per_time_step.smart_kinetics_estimate;
    j["computing_costs_per_time_step"]["smart_kinetics_search"] = obj.computing_costs_per_time_step.smart_kinetics_search;
    j["computing_costs_per_time_step"]["smart_kinetics_error_control"] = obj.computing_costs_per_time_step.smart_kinetics_error_control;
    j["computing_costs_per_time_step"]["smart_kinetics_taylor"] = obj.computing_costs_per_time_step.smart_kinetics_taylor;

    j["computing_costs_per_time_step"]["smart_kinetics_learn"] = obj.computing_costs_per_time_step.smart_kinetics_learn;
    j["computing_costs_per_time_step"]["smart_kinetics_chemical_properties"] = obj.computing_costs_per_time_step.smart_kinetics_chemical_properties;
    j["computing_costs_per_time_step"]["smart_kinetics_equilibration"] = obj.computing_costs_per_time_step.smart_kinetics_equilibration;

    j["computing_costs_per_time_step"]["smart_equilibrium"] = obj.computing_costs_per_time_step.smart_equilibrium;
    j["computing_costs_per_time_step"]["smart_equilibrium_with_ideal_search"] = obj.computing_costs_per_time_step.smart_equilibrium_with_ideal_search;

    j["computing_costs_per_time_step"]["smart_equilibrium_estimate"] = obj.computing_costs_per_time_step.smart_equilibrium_estimate;
    j["computing_costs_per_time_step"]["smart_equilibrium_search"] = obj.computing_costs_per_time_step.smart_equilibrium_search;
    j["computing_costs_per_time_step"]["smart_equilibrium_error_control"] = obj.computing_costs_per_time_step.smart_equilibrium_error_control;
    j["computing_costs_per_time_step"]["smart_equilibrium_taylor"] = obj.computing_costs_per_time_step.smart_equilibrium_taylor;
    j["computing_costs_per_time_step"]["smart_equilibrium_database_priority_update"] = obj.computing_costs_per_time_step.smart_equilibrium_database_priority_update;

    j["computing_costs_per_time_step"]["smart_equilibrium_learn"] = obj.computing_costs_per_time_step.smart_equilibrium_learn;
    j["computing_costs_per_time_step"]["smart_equilibrium_gibbs_energy_minimization"] = obj.computing_costs_per_time_step.smart_equilibrium_gibbs_energy_minimization;
    j["computing_costs_per_time_step"]["smart_equilibrium_chemical_properties"] = obj.computing_costs_per_time_step.smart_equilibrium_chemical_properties;
    j["computing_costs_per_time_step"]["smart_equilibrium_sensitivity_matrix"] = obj.computing_costs_per_time_step.smart_equilibrium_sensitivity_matrix;
    j["computing_costs_per_time_step"]["smart_equilibrium_error_control_matrices"] = obj.computing_costs_per_time_step.smart_equilibrium_error_control_matrices;
    j["computing_costs_per_time_step"]["smart_equilibrium_storage"] = obj.computing_costs_per_time_step.smart_equilibrium_storage;

}

void from_json(const json& j, ReactiveTransportAnalysis& obj) {
    j.at("transport").at("timing").get_to(obj.transport.timing);

    j.at("equilibrium").at("timing").get_to(obj.equilibrium.timing);

    j.at("kinetics").at("timing").get_to(obj.kinetics.timing);

    j.at("smart_equilibrium").at("timing").get_to(obj.smart_equilibrium.timing);
    j.at("num_equilibrium_calculations").get_to(obj.smart_equilibrium.num_equilibrium_calculations);
    j.at("num_smart_equilibrium_accepted_estimates").get_to(obj.smart_equilibrium.num_smart_equilibrium_accepted_estimates);
    j.at("num_smart_equilibrium_required_learnings").get_to(obj.smart_equilibrium.num_smart_equilibrium_required_learnings);
    j.at("smart_equilibrium_estimate_acceptance_rate").get_to(obj.smart_equilibrium.smart_equilibrium_estimate_acceptance_rate);
    j.at("smart_equilibrium_cells_where_learning_was_required_at_step").get_to(obj.smart_equilibrium.cells_where_learning_was_required_at_step);

    j.at("smart_kinetics").at("timing").get_to(obj.smart_kinetics.timing);
    j.at("num_kinetics_calculations").get_to(obj.smart_kinetics.num_kinetics_calculations);
    j.at("num_smart_kinetics_accepted_estimates").get_to(obj.smart_kinetics.num_smart_kinetics_accepted_estimates);
    j.at("num_smart_kinetics_required_learnings").get_to(obj.smart_kinetics.num_smart_kinetics_required_learnings);
    j.at("smart_kinetics_estimate_acceptance_rate").get_to(obj.smart_kinetics.smart_kinetics_estimate_acceptance_rate);
    j.at("smart_kinetics_cells_where_learning_was_required_at_step").get_to(obj.smart_kinetics.cells_where_learning_was_required_at_step);

    j.at("computing_costs_per_time_step").at("transport").get_to(obj.computing_costs_per_time_step.transport);
    j.at("computing_costs_per_time_step").at("equilibrium").get_to(obj.computing_costs_per_time_step.equilibrium);
    j.at("computing_costs_per_time_step").at("kinetics").get_to(obj.computing_costs_per_time_step.kinetics);

    j.at("computing_costs_per_time_step").at("smart_kinetics").get_to(obj.computing_costs_per_time_step.smart_kinetics);
    j.at("computing_costs_per_time_step").at("smart_kinetics_estimate").get_to(obj.computing_costs_per_time_step.smart_kinetics_estimate);
    j.at("computing_costs_per_time_step").at("smart_kinetics_search").get_to(obj.computing_costs_per_time_step.smart_kinetics_search);
    j.at("computing_costs_per_time_step").at("smart_kinetics_error_taylor").get_to(obj.computing_costs_per_time_step.smart_kinetics_taylor);
    j.at("computing_costs_per_time_step").at("smart_kinetics_error_control").get_to(obj.computing_costs_per_time_step.smart_kinetics_error_control);
    j.at("computing_costs_per_time_step").at("smart_kinetics_learn").get_to(obj.computing_costs_per_time_step.smart_kinetics_learn);
    j.at("computing_costs_per_time_step").at("smart_kinetics_chemical_properties").get_to(obj.computing_costs_per_time_step.smart_kinetics_chemical_properties);
    j.at("computing_costs_per_time_step").at("smart_kinetics_equilibration").get_to(obj.computing_costs_per_time_step.smart_kinetics_equilibration);

    j.at("computing_costs_per_time_step").at("smart_equilibrium").get_to(obj.computing_costs_per_time_step.smart_equilibrium);
    j.at("computing_costs_per_time_step").at("smart_equilibrium_with_ideal_search").get_to(obj.computing_costs_per_time_step.smart_equilibrium_with_ideal_search);

    j.at("computing_costs_per_time_step").at("smart_equilibrium_estimate").get_to(obj.computing_costs_per_time_step.smart_equilibrium_estimate);
    j.at("computing_costs_per_time_step").at("smart_equilibrium_search").get_to(obj.computing_costs_per_time_step.smart_equilibrium_search);
    j.at("computing_costs_per_time_step").at("smart_equilibrium_error_control").get_to(obj.computing_costs_per_time_step.smart_equilibrium_error_control);
    j.at("computing_costs_per_time_step").at("smart_equilibrium_error_taylor").get_to(obj.computing_costs_per_time_step.smart_equilibrium_taylor);
    j.at("computing_costs_per_time_step").at("smart_equilibrium_error_database_priority_update").get_to(obj.computing_costs_per_time_step.smart_equilibrium_database_priority_update);

    j.at("computing_costs_per_time_step").at("smart_equilibrium_learn").get_to(obj.computing_costs_per_time_step.smart_equilibrium_learn);
    j.at("computing_costs_per_time_step").at("smart_equilibrium_gibbs_energy_minimization").get_to(obj.computing_costs_per_time_step.smart_equilibrium_gibbs_energy_minimization);
    j.at("computing_costs_per_time_step").at("smart_equilibrium_chemical_properties").get_to(obj.computing_costs_per_time_step.smart_equilibrium_chemical_properties);
    j.at("computing_costs_per_time_step").at("smart_equilibrium_sensitivity_matrix").get_to(obj.computing_costs_per_time_step.smart_equilibrium_sensitivity_matrix);
    j.at("computing_costs_per_time_step").at("smart_equilibrium_error_control_matrices").get_to(obj.computing_costs_per_time_step.smart_equilibrium_error_control_matrices);
    j.at("computing_costs_per_time_step").at("smart_equilibrium_storage").get_to(obj.computing_costs_per_time_step.smart_equilibrium_storage);
}

} // namespace Reaktoro
