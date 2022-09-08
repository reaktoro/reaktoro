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

#include "ReactionRateModels.hpp"

// Reaktoro includes
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Models/ReactionRateModels.hpp>

namespace Reaktoro {

//======================================================================
// ReactionRateModelParams Types
//======================================================================

REAKTORO_DATA_ENCODE_DECLARE(ReactionRateModelParamsPalandriKharaka::Mechanism);
REAKTORO_DATA_DECODE_DECLARE(ReactionRateModelParamsPalandriKharaka::Mechanism);

REAKTORO_DATA_ENCODE_DEFINE(ReactionRateModelParamsPalandriKharaka::Mechanism)
{
    data["lgk"] = obj.lgk;
    data["E"] = obj.E;
    data["p"] = obj.p;
    data["q"] = obj.q;
    for(auto catalyst : obj.catalysts)
    {
        const auto key = catalyst.property + "(" + catalyst.formula + ")";
        data[key] = catalyst.power;
    }
}

REAKTORO_DATA_DECODE_DEFINE(ReactionRateModelParamsPalandriKharaka::Mechanism)
{
    data.required("lgk").to(obj.lgk);
    data.required("E").to(obj.E);
    data.optional("p").to(obj.p);
    data.optional("q").to(obj.q);

    // Collect all catalyst properties and their power
    for(auto const& [key, value] : data.asDict())
    {
        if(!oneof(key[0], 'a', 'P'))
            continue;
        errorif(key.size() <= 3, "Expecting a chemical formula inside `a()` or `P()`, such as `a(H+)`, `P(CO2)`.");
        errorif(key[1] != '(' || key.back() != ')', "Expecting ( and ) as in `a(H+)`, `a(Fe+3)`, `P(CO2)`.");
        const auto formula = key.substr(2, key.size() - 3); // exclude first two chars and last
        const auto property = key.substr(0, 1);
        const auto power = value.asFloat();
        obj.catalysts.push_back({ formula, property, power });
    }
}

REAKTORO_DATA_ENCODE_DEFINE(ReactionRateModelParamsPalandriKharaka)
{
    data["Mineral"] = join(obj.names, " ");
    for(auto mechanism : obj.mechanisms)
        data["Mechanisms"][mechanism.name] = mechanism;
}

REAKTORO_DATA_DECODE_DEFINE(ReactionRateModelParamsPalandriKharaka)
{
    obj.names = split(data.required("Mineral").asString(), " ");
    for(auto const& [name, mechanism] : data.required("Mechanisms").asDict())
    {
        obj.mechanisms.push_back(mechanism.as<ReactionRateModelParamsPalandriKharaka::Mechanism>());
        obj.mechanisms.back().name = name;
    }
}

} // namespace Reaktoro
