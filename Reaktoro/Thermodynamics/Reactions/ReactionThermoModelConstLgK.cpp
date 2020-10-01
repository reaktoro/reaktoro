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

#include "ReactionThermoModelConstLgK.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>

namespace Reaktoro {

auto ReactionThermoModelConstLgK(real lgK0) -> ReactionThermoPropsFn
{
    return [=](real T, real P) -> ReactionThermoProps
    {
        const auto R = universalGasConstant;
        const auto lnK0 = lgK0 * ln10;
        ReactionThermoProps props;
        props.dG0 = -R*T*lnK0;
        props.dH0 = 0.0;
        return props;
    };
}

} // namespace Reaktoro
