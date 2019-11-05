// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/AqueousMixture.hpp>

namespace Reaktoro {

auto aqueousActivityModelDrummondCO2(const AqueousMixture& mixture) -> AqueousActivityModel
{
    AqueousActivityModel f = [=](const AqueousMixtureState& state) {
        // Calculate the activity coefficient of CO2(aq)
        const ThermoScalar& T = state.T;
        const ThermoScalar c1 = -1.0312 + 1.2806e-3 * T + 255.9 / T;
        const ThermoScalar c2 = 0.4445 - 1.6060e-3 * T;

        // The stoichiometric ionic strength of the aqueous mixture
        const auto& I = state.Is;

        // The ln activity coefficient of CO2(aq)
        ChemicalScalar ln_gCO2 = c1 * I - c2 * I / (I + 1);

        return ln_gCO2;
    };

    return f;
}

} // namespace Reaktoro
