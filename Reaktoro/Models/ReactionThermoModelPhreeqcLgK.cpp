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

#include "ReactionThermoModelPhreeqcLgK.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Serialization/Serialization.YAML.hpp>

namespace Reaktoro {

/// Return a Params object containing all Param objects in @p params.
auto extractParams(const ReactionThermoModelParamsPhreeqcLgK& params) -> Params
{
    const auto& [A1, A2, A3, A4, A5, A6, Pr] = params;
    return {A1, A2, A3, A4, A5, A6};
}

/// Return a ModelSerializer for given model parameters in @p params.
auto createModelSerializer(const ReactionThermoModelParamsPhreeqcLgK& params) -> ModelSerializer
{
    return [=]()
    {
        yaml node;
        node["PhreeqcLgK"] = params;
        return node;
    };
}

auto ReactionThermoModelPhreeqcLgK(const ReactionThermoModelParamsPhreeqcLgK& params) -> ReactionThermoModel
{
    auto evalfn = [=](ReactionThermoProps& props, ReactionThermoArgs args)
    {
        // Unpack the arguments for the evaluation of this model
        const auto& [T, P, dV0] = args;

        // Unpack the model parameters
        const auto& [A1, A2, A3, A4, A5, A6, Pr] = params;

        const auto R = universalGasConstant;
        const auto T2 = T*T;
        const auto T3 = T*T2;
        const auto dE = dV0 * (P - Pr); // delta energy (in J/mol)
        props.dG0 = -R*T * (A1 + A2*T + A3/T + A4*log10(T) + A5/T2 + A6*T2)*ln10 + dE;
        props.dH0 = R * (A2*T2 - A3 + A4*T/ln10 - 2*A5/T + 2*A6*T3)*ln10 + dE;
    };

    return ReactionThermoModel(evalfn, extractParams(params), createModelSerializer(params));
}

} // namespace Reaktoro
