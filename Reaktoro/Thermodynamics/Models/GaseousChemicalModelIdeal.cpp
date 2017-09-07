// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "GaseousChemicalModelIdeal.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/GaseousMixture.hpp>

namespace Reaktoro {

auto gaseousChemicalModelIdeal(const GaseousMixture& mixture) -> PhaseChemicalModel
{
    // The state of the gaseous mixture
    GaseousMixtureState state;

    // Define the chemical model function of the gaseous phase
    PhaseChemicalModel model = [=](PhaseChemicalModelResult& res, Temperature T, Pressure P, const Vector& n) mutable
    {
        // Evaluate the state of the gaseous mixture
        state = mixture.state(T, P, n);

        // Calculate pressure in bar
        const ThermoScalar Pbar = 1e-5 * Pressure(P);

        // The ln of pressure in units of bar
        const ThermoScalar ln_Pbar = log(Pbar);

        // The result of the ideal model
        res.ln_activities = log(state.x) + ln_Pbar;
    };

    return model;
}

} // namespace Reaktoro



