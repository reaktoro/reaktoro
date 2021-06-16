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

#include "ReactionThermoModelGemsLgK.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Serialization/Models.YAML.hpp>

namespace Reaktoro {

/// Return a Params object containing all Param objects in @p params.
auto extractParams(const ReactionThermoModelParamsGemsLgK& params) -> Params
{
    const auto& [A0, A1, A2, A3, A4, A5, A6, Pr] = params;
    return {A0, A1, A2, A3, A4, A5, A6};
}

/// Return a ModelSerializer for given model parameters in @p params.
auto createModelSerializer(const ReactionThermoModelParamsGemsLgK& params) -> ModelSerializer
{
    return [=]()
    {
        yaml node;
        node["GemsLgK"] = params;
        return node;
    };
}

auto ReactionThermoModelGemsLgK(const ReactionThermoModelParamsGemsLgK& params) -> ReactionThermoModel
{
    auto evalfn = [=](ReactionThermoProps& props, ReactionThermoArgs args)
    {
        // Unpack the arguments for the evaluation of this model
        const auto& [T, P, dV0] = args;

        // Unpack the model parameters
        const auto& [A0, A1, A2, A3, A4, A5, A6, Pr] = params;

        const auto R = universalGasConstant;
        const auto T2 = T*T;
        const auto T3 = T*T2;
        const auto T05 = sqrt(T);
        const auto dE = dV0 * (P - Pr); // delta energy (in J/mol)

        props.dG0 = -R*T * (A0 + A1*T + A2/T + A3*log(T) + A4/T2 + A5*T2 + A6/T05)*ln10 + dE;
        props.dH0 = R * (A1*T2 - A2 + A3*T - 2*A4/T + 2*A5*T3 - 0.5*A6*T05)*ln10 + dE;
    };

    return ReactionThermoModel(evalfn, extractParams(params), createModelSerializer(params));
}

} // namespace Reaktoro
