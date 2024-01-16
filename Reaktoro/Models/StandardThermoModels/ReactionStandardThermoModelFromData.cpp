// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

#include "ReactionStandardThermoModelFromData.hpp"

// Reaktoro includes
#include <Reaktoro/Models/StandardThermoModels.hpp>
#include <Reaktoro/Serialization/Models/StandardThermoModels.hpp>

namespace Reaktoro {

auto ReactionStandardThermoModelFromData(Data const& data) -> ReactionStandardThermoModel
{
    errorif(!data.isDict(),
        "Expecting a dictionary containing a single key-value pair in the Data object:\n", data.repr());

    const auto model = data.asDict();

    errorif(model.size() != 1,
        "Expecting only one key-value pair in the Data object:\n", data.repr());

    auto const& name = model.front().first;
    auto const& params = model.front().second;

    if(name == "ConstLgK")
        return ReactionStandardThermoModelConstLgK(params.as<ReactionStandardThermoModelParamsConstLgK>());
    if(name == "GemsLgK")
        return ReactionStandardThermoModelGemsLgK(params.as<ReactionStandardThermoModelParamsGemsLgK>());
    if(name == "PhreeqcLgK")
        return ReactionStandardThermoModelPhreeqcLgK(params.as<ReactionStandardThermoModelParamsPhreeqcLgK>());
    if(name == "VantHoff")
        return ReactionStandardThermoModelVantHoff(params.as<ReactionStandardThermoModelParamsVantHoff>());

    errorif(true, "Cannot create a ReactionStandardThermoModel object with "
        "unsupported model name `", name, "` in Data object:\n", data.repr());

    return {};
}

} // namespace Reaktoro
