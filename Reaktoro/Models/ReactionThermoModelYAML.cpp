// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

#include "ReactionThermoModelYAML.hpp"

// Reaktoro includes
#include <Reaktoro/Models/ReactionThermoModelConstLgK.hpp>
#include <Reaktoro/Models/ReactionThermoModelGemsLgK.hpp>
#include <Reaktoro/Models/ReactionThermoModelPhreeqcLgK.hpp>
#include <Reaktoro/Models/ReactionThermoModelVantHoff.hpp>
#include <Reaktoro/Serialization/Models.YAML.hpp>

namespace Reaktoro {

auto ReactionThermoModelYAML(const yaml& node) -> ReactionThermoModel
{
    errorif(!(node.IsMap() && node.size() == 1),
        "Expecting the following yaml node to contain a single "
        "StandardThermoModel child:\n", node.repr());

    auto model = node.begin()->first.as<String>();
    yaml params = node.begin()->second;

    if(model == "ConstLgK")
        return ReactionThermoModelConstLgK(params);
    if(model == "GemsLgK")
        return ReactionThermoModelGemsLgK(params);
    if(model == "PhreeqcLgK")
        return ReactionThermoModelPhreeqcLgK(params);
    if(model == "VantHoff")
        return ReactionThermoModelVantHoff(params);
    errorif(true, "Cannot create a ReactionThermoModel with "
        "unsupported model name `", model, "` in yaml node:\n", node.repr());
}

} // namespace Reaktoro
