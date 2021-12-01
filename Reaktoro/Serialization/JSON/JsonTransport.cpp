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
#include <Reaktoro/Serialization/JSON/JsonEquilibrium.hpp>
#include <Reaktoro/Transport/TransportResult.hpp>
#include <Reaktoro/Transport/ReactiveTransportAnalysis.hpp>

namespace Reaktoro {

void to_json(json& j, const TransportTiming& obj)
{
    j["step"] = obj.step;
    j["matrix_equation_assembly"] = obj.matrix_equation_assembly;
    j["matrix_equation_solve"] = obj.matrix_equation_solve;
}

void from_json(const json& j, TransportTiming& obj)
{
    j.at("step").get_to(obj.step);
    j.at("matrix_equation_assembly").get_to(obj.matrix_equation_assembly);
    j.at("matrix_equation_solve").get_to(obj.matrix_equation_solve);
}

void to_json(json& j, const TransportResult& obj)
{
    j["timing"] = obj.timing;
}

void from_json(const json& j, TransportResult& obj)
{
    j.at("timing").get_to(obj.timing);
}

void to_json(json& j, const ReactiveTransportAnalysis& obj)
{
    j["transport"] = json::object();
    j["transport"]["timing"] = obj.transport.timing;

    j["equilibrium"] = json::object();
    j["equilibrium"]["timing"] = obj.equilibrium.timing;

    j["computing_costs_per_time_step"] = json::object();
    j["computing_costs_per_time_step"]["transport"] = obj.computing_costs_per_time_step.transport;
    j["computing_costs_per_time_step"]["equilibrium"] = obj.computing_costs_per_time_step.equilibrium;
}

void from_json(const json& j, ReactiveTransportAnalysis& obj)
{
    j.at("transport").at("timing").get_to(obj.transport.timing);

    j.at("equilibrium").at("timing").get_to(obj.equilibrium.timing);

    j.at("computing_costs_per_time_step").at("transport").get_to(obj.computing_costs_per_time_step.transport);
    j.at("computing_costs_per_time_step").at("equilibrium").get_to(obj.computing_costs_per_time_step.equilibrium);

}

} // namespace Reaktoro
