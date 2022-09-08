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

#include "ReactionStandardThermoModelYAML.hpp"

// Reaktoro includes
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelConstLgK.hpp>
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelGemsLgK.hpp>
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelPhreeqcLgK.hpp>
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelVantHoff.hpp>
#include <Reaktoro/Serialization/Models.YAML.hpp>

namespace Reaktoro {

auto ReactionStandardThermoModelYAML(const yaml& node) -> ReactionStandardThermoModel
{
    errorif(!(node.IsMap() && node.size() == 1),
        "Expecting the following yaml node to contain a single "
        "StandardThermoModel child:\n", node.repr());

    auto model = node.begin()->first.as<String>();
    yaml params = node.begin()->second;

    if(model == "ConstLgK")
        return ReactionStandardThermoModelConstLgK(params);
    if(model == "GemsLgK")
        return ReactionStandardThermoModelGemsLgK(params);
    if(model == "PhreeqcLgK")
        return ReactionStandardThermoModelPhreeqcLgK(params);
    if(model == "VantHoff")
        return ReactionStandardThermoModelVantHoff(params);
    errorif(true, "Cannot create a ReactionStandardThermoModel with "
        "unsupported model name `", model, "` in yaml node:\n", node.repr());

    return {};
}

} // namespace Reaktoro
