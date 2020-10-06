// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

#include "ReactionThermoModelAnalyticalPHREEQC.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>

namespace Reaktoro {

auto ReactionThermoModelAnalyticalPHREEQC(real A1, real A2, real A3, real A4, real A5, real A6) -> ReactionThermoModel
{
    auto creatorfn = [](const Params& params)
    {
        const real A1 = params.get("A1");
        const real A2 = params.get("A2");
        const real A3 = params.get("A3");
        const real A4 = params.get("A4");
        const real A5 = params.get("A5");
        const real A6 = params.get("A6");

        return [=](ReactionThermoProps& props, ReactionThermoArgs args)
        {
            ReactionThermoArgsDecl(args);
            const auto R = universalGasConstant;
            const auto T2 = T*T;
            const auto T3 = T*T2;
            props.dG0 = -R*T * (A1 + A2*T + A3/T + A4*log10(T) + A5/T2 + A6*T2) * ln10;
            props.dH0 = R * (A2*T2 - A3 + A4*T/ln10 - 2*A5/T + 2*A6*T3) * ln10;
        };
    };

    Params params;
    params.set("A1", A1);
    params.set("A2", A2);
    params.set("A3", A3);
    params.set("A4", A4);
    params.set("A5", A5);
    params.set("A6", A6);

    return ReactionThermoModel(creatorfn, params);
}

} // namespace Reaktoro
