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

#include "ReactionThermoModelPressureCorrection.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>

namespace Reaktoro {

auto ReactionThermoModelPressureCorrection(Param Pr) -> ReactionThermoModel
{
    auto evalfn = [=](ReactionThermoProps& props, ReactionThermoArgs args)
    {
        // Unpack the arguments for the evaluation of this model
        const auto& [T, P, dV0] = args;

        const auto dE = (P - Pr) * dV0; // delta energy (in J/mol)
        props.dG0 += dE;
        props.dH0 += dE;
    };

    Vec<Param> params = { Pr };

    return ReactionThermoModel(evalfn, params);
}

} // namespace Reaktoro
