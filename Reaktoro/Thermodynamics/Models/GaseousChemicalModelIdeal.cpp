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

#include "GaseousChemicalModelIdeal.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/GaseousMixture.hpp>

namespace Reaktoro {

auto gaseousChemicalModelIdeal(const GaseousMixture& mixture) -> PhaseChemicalModel
{
    const unsigned nspecies = mixture.numSpecies();

    PhaseChemicalModel f = [=](double T, double P, const Vector& n)
    {
        // Calculate the state of the mixture
        GaseousMixtureState state = mixture.state(T, P, n);

        // Calculate pressure in bar
        const ThermoScalar Pbar = 1e-5 * Pressure(P);

        // The ln of pressure in units of bar
        const ThermoScalar ln_Pbar = log(Pbar);

        // The result of the ideal model
        PhaseChemicalModelResult res(nspecies);
        res.ln_activity_constants = ln_Pbar;
        res.ln_activities = log(state.x) + ln_Pbar;

        return res;
    };

    return f;
}

} // namespace Reaktoro



