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
    const Index iH2O = mixture.indexWater();

    // The state of the aqueous mixture
    AqueousMixtureState state;

    PhaseChemicalModel f = [=](PhaseChemicalModelResult& res, Temperature T, Pressure P, VectorConstRef n) mutable
    {
        // Evaluate the state of the aqueous mixture
        state = mixture.state(T, P, n);

        // The ln of water mole fraction
        ChemicalScalar ln_xw = log(state.x[iH2O]);

        // Set the activity coefficients of the aqueous species
        res.ln_activity_coefficients = ln_xw;
        res.ln_activity_coefficients[iH2O] = 0.0;

        // Set the activities of the aqueous species
        res.ln_activities = res.ln_activity_coefficients + log(state.m);
        res.ln_activities[iH2O] = ln_xw;
    };

    return f;
}

} // namespace Reaktoro



