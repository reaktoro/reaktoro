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

#include "ReactionStandardThermoModelPhreeqcLgK.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Serialization/Models/StandardThermoModels.hpp>

namespace Reaktoro {

auto ReactionStandardThermoModelPhreeqcLgK(const ReactionStandardThermoModelParamsPhreeqcLgK& params) -> ReactionStandardThermoModel
{
    auto evalfn = [=](ReactionStandardThermoProps& props, ReactionStandardThermoModelArgs args)
    {
        // Unpack the arguments for the evaluation of this model
        const auto& [T, P, dV0] = args;

        // Unpack the model parameters
        const auto& [A1, A2, A3, A4, A5, A6, Pr] = params;

        const auto dVT0 = 0.0; // TODO: Consider dVT0 in ReactionStandardThermoModelArgs

        const auto R = universalGasConstant;
        const auto T2 = T*T;
        const auto T3 = T*T2;
        const auto dE = dV0 * (P - Pr); // delta energy (in J/mol)
        const auto dET = dVT0 * (P - Pr);
        props.dG0  = -R*T * (A1 + A2*T + A3/T + A4*log10(T) + A5/T2 + A6*T2)*ln10 + dE;
        props.dH0  = R * (A2*T2 - A3 + A4*T/ln10 - 2*A5/T + 2*A6*T3)*ln10 + dE;
        props.dCp0 = R * (2*A2*T + A4/ln10 + 2*A5/T2 + 6*A6*T2)*ln10 + dET;
    };

    Data paramsdata;
    paramsdata["PhreeqcLgK"] = params;

    return ReactionStandardThermoModel(evalfn, paramsdata);
}

} // namespace Reaktoro
