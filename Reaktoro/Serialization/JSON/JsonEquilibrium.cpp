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

#include "JsonEquilibrium.hpp"

// Reaktoro includes
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>

namespace Reaktoro {

void to_json(json& j, const EquilibriumResult& obj) {
    j["timing"] = obj.timing;
}

void from_json(const json& j, EquilibriumResult& obj) {
    j.at("timing").get_to(obj.timing);
}

void to_json(json& j, const EquilibriumTiming& obj) {
    j["solve"] = obj.solve;
    j["standard_thermodynamic_properties"] = obj.standard_thermodynamic_properties;
    j["chemical_properties"] = obj.chemical_properties;
}

void from_json(const json& j, EquilibriumTiming& obj) {
    j.at("solve").get_to(obj.solve);
    j.at("standard_thermodynamic_properties").get_to(obj.standard_thermodynamic_properties);
    j.at("chemical_properties").get_to(obj.chemical_properties);
}

} // namespace Reaktoro
