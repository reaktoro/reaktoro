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

#include "AqueousActivityModelDrummondCO2.hpp"

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Mixtures/AqueousMixture.hpp>

namespace Reaktoro {

auto aqueousActivityModelDrummondCO2(const AqueousMixture& mixture) -> AqueousActivityModel
{
    AqueousActivityModel f = [=](const AqueousMixtureState& state)
    {
        const auto T = state.T;
        const auto c1 = -1.0312 + 1.2806e-3*T + 255.9/T;
        const auto c2 =  0.4445 - 1.6060e-3*T;
        const auto I = state.Is;
        return c1 * I - c2 * I/(I + 1);
    };

    return f;
}

} // namespace Reaktoro
