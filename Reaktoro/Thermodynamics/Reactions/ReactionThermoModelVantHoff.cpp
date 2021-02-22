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

#include "ReactionThermoModelVantHoff.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>

namespace Reaktoro {

auto ReactionThermoModelVantHoff(Param lgK0, Param dH0, Param Tr) -> ReactionThermoModel
{
    auto evalfn = [=](ReactionThermoProps& props, ReactionThermoArgs args)
    {
        ReactionThermoArgsDecl(args);
        const auto R = universalGasConstant;
        const auto RT = R*T;
        const auto lnK0 = lgK0 * ln10;
        const auto lnK = lnK0 - dH0 * (Tr - T)/(RT*Tr);
        props.dG0 = -RT * lnK;
        props.dH0 = dH0;
    };

    Params params = { lgK0, dH0 };

    return ReactionThermoModel(evalfn, params);
}

} // namespace Reaktoro
