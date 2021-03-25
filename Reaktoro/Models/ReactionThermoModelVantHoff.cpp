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

#include "ReactionThermoModelVantHoff.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Serialization/SerializationYAML.hpp>

namespace Reaktoro {

/// Return a Params object containing all Param objects in @p params.
auto extractParams(const ReactionThermoModelParamsVantHoff& params) -> Params
{
    const auto& [lgKr, dHr, Tr, Pr] = params;
    return {lgKr, dHr};
}

/// Return a ModelSerializer for given model parameters in @p params.
auto createModelSerializer(const ReactionThermoModelParamsVantHoff& params) -> ModelSerializer
{
    return [=]()
    {
        yaml node;
        node["VantHoff"] = params;
        return node;
    };
}

auto ReactionThermoModelVantHoff(const ReactionThermoModelParamsVantHoff& params) -> ReactionThermoModel
{
    auto evalfn = [=](ReactionThermoProps& props, ReactionThermoArgs args)
    {
        // Unpack the arguments for the evaluation of this model
        const auto& [T, P, dV0] = args;

        // Unpack the model parameters
        const auto& [lgKr, dHr, Tr, Pr] = params;

        const auto R = universalGasConstant;
        const auto RT = R*T;
        const auto lnKr = lgKr * ln10;
        const auto lnK = lnKr - dHr * (Tr - T)/(RT*Tr);
        const auto dE = dV0 * (P - Pr); // delta energy (in J/mol)
        props.dG0 = -RT * lnK + dE;
        props.dH0 = dHr + dE;
    };

    return ReactionThermoModel(evalfn, extractParams(params), createModelSerializer(params));
}

} // namespace Reaktoro
