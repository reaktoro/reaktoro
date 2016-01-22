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

#include "AqueousChemicalModelIdeal.hpp"

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Mixtures/AqueousMixture.hpp>

namespace Reaktoro {

auto aqueousChemicalModelIdeal(const AqueousMixture& mixture) -> PhaseChemicalModel
{
    const unsigned nspecies = mixture.numSpecies();
    const Index iH2O = mixture.indexWater();

    PhaseChemicalModel f = [=](double T, double P, const Vector& n)
    {
        // Calculate state of the mixture
        const AqueousMixtureState state = mixture.state(T, P, n);

        PhaseChemicalModelResult res(nspecies);
        res.ln_activities = log(state.m);
        res.ln_activities[iH2O] = log(state.x[iH2O]);
        return res;
    };

    return f;
}

} // namespace Reaktoro



