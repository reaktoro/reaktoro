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

#include "AqueousChemicalModelIdeal.hpp"

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Mixtures/AqueousMixture.hpp>

namespace Reaktoro {

auto aqueousChemicalModelIdeal(const AqueousMixture& mixture)-> ActivityModelFn
{
    const Index iH2O = mixture.indexWater();

    // The state of the aqueous mixture
    AqueousMixtureState state;

    ActivityModelFn f = [=](ActivityProps res, real T, real P, ArrayXrConstRef x) mutable
    {
        using std::log;

        // Evaluate the state of the aqueous mixture
        state = mixture.state(T, P, x);

        // The ln of water mole fraction
        real ln_xw = log(x[iH2O]);

        // Set the activity coefficients of the aqueous species
        res.ln_g = ln_xw;
        res.ln_g[iH2O] = 0.0;

        // Set the activities of the aqueous species
        res.ln_a = res.ln_g + state.m.log();
        res.ln_a[iH2O] = ln_xw;
    };

    return f;
}

} // namespace Reaktoro



