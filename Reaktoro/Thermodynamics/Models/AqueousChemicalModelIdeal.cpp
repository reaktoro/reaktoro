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

        // The ln of water molar fraction
        ChemicalScalar ln_xw = log(state.x[iH2O]);

        // The result of the phase chemical model
        PhaseChemicalModelResult res(nspecies);

        // Set the activity coefficients of the aqueous species
        res.ln_activity_coefficients = ln_xw;
        res.ln_activity_coefficients[iH2O] = 0.0;

        // Set the activities of the aqueous species
        res.ln_activities = res.ln_activity_coefficients + log(state.m);
        res.ln_activities[iH2O] = ln_xw;

        // Set the activity constants of the aqueous species
        res.ln_activity_constants = std::log(55.508472);
        res.ln_activity_constants[iH2O] = 0.0;

        return res;
    };

    return f;
}

} // namespace Reaktoro



