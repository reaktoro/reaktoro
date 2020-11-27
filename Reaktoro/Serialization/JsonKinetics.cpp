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

#include "JsonKinetics.hpp"

// Reaktoro includes
#include <Reaktoro/Kinetics/KineticResult.hpp>
#include <Reaktoro/Kinetics/SmartKineticResult.hpp>

namespace Reaktoro {

void to_json(json& j, const KineticResult& obj) {
    j["timing"] = obj.timing;
    // j["optimum"] = obj.optimum;
}

void from_json(const json& j, KineticResult& obj) {
    j.at("timing").get_to(obj.timing);
    // j.at("optimum").get_to(obj.optimum);
}

void to_json(json& j, const KineticTiming& obj) {
    j["solve"] = obj.solve;
    j["initialize"] = obj.initialize;
    j["integrate"] = obj.integrate;
    j["chemical_properties"] = obj.integrate_chemical_properties;
    j["integrate_reaction_rates"] = obj.integrate_reaction_rates;
    j["equilibration"] = obj.integrate_equilibration;

}

void from_json(const json& j, KineticTiming& obj) {
    j.at("solve").get_to(obj.solve);
    j.at("initialize").get_to(obj.initialize);
    j.at("integrate").get_to(obj.integrate);
    j.at("chemical_properties").get_to(obj.integrate_chemical_properties);
    j.at("integrate_reaction_rates").get_to(obj.integrate_reaction_rates);
    j.at("equilibration").get_to(obj.integrate_equilibration);
}

void to_json(json& j, const SmartKineticResult& obj) {
    j["estimate"] = json::object();
    j["estimate"]["accepted"] = obj.estimate.accepted;
    j["estimate"]["failed_with_species"] = obj.estimate.failed_with_species;
    j["estimate"]["failed_with_amount"] = obj.estimate.failed_with_amount;
    j["estimate"]["failed_with_chemical_potential"] = obj.estimate.failed_with_chemical_potential;
    //j["learn"] = json::object();
    //j["learn"]["gibbs_energy_minimization"] = obj.learning.gibbs_energy_minimization;
    j["timing"] = obj.timing;
}

void from_json(const json& j, SmartKineticResult& obj) {
    j.at("estimate").at("accepted").get_to(obj.estimate.accepted);
    j.at("estimate").at("failed_with_species").get_to(obj.estimate.failed_with_species);
    j.at("estimate").at("failed_with_amount").get_to(obj.estimate.failed_with_amount);
    j.at("estimate").at("failed_with_chemical_potential").get_to(obj.estimate.failed_with_chemical_potential);
    //j.at("learn").at("gibbs_energy_minimization").get_to(obj.learning.gibbs_energy_minimization);
    j.at("timing").get_to(obj.timing);
}

void to_json(json& j, const SmartKineticTiming& obj) {
    j["initialize"] = obj.initialize;
    j["solve"] = obj.solve;
    j["learn"] = obj.learn;
    j["learn_integrate"] = obj.learn_integration;
    j["learn_chemical_properties"] = obj.learn_chemical_properties;
    j["learn_reaction_rates"] = obj.learn_reaction_rates;
    j["learn_equilibration"] = obj.learn_equilibration;
    j["estimate"] = obj.estimate;
    j["estimate_search"] = obj.estimate_search;
    j["estimate_taylor"] = obj.estimate_taylor;
    j["estimate_error_control"] = obj.estimate_error_control;
}

void from_json(const json& j, SmartKineticTiming& obj) {
    j.at("initialize").get_to(obj.initialize);
    j.at("solve").get_to(obj.solve);
    j.at("learn").get_to(obj.learn);
    j.at("learn_integrate").get_to(obj.learn_integration);
    j.at("learn_chemical_properties").get_to(obj.learn_chemical_properties);
    j.at("learn_reaction_rates").get_to(obj.learn_reaction_rates);
    j.at("learn_equilibration").get_to(obj.learn_equilibration);
    j.at("estimate").get_to(obj.estimate);
    j.at("estimate_search").get_to(obj.estimate_search);
    j.at("estimate_taylor").get_to(obj.estimate_taylor);
    j.at("estimate_error_control").get_to(obj.estimate_error_control);
}

} // namespace Reaktoro
