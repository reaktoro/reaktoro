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

#include "AqueousActivityModelRumpfCO2.hpp"

// Reaktoro includes
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/AqueousMixture.hpp>

namespace Reaktoro {

auto aqueousActivityModelRumpfCO2(const AqueousMixture& mixture) -> AqueousActivityModel
{
    // The number of speciesn and charged species
    const unsigned nspecies = mixture.numSpecies();
    const unsigned nions = mixture.numChargedSpecies();

    // The local indices of some charged species among all charged species
    const Index iNa  = mixture.indexChargedSpeciesAny(alternativeChargedSpeciesNames("Na+"));  // Na+, Na[+]
    const Index iK   = mixture.indexChargedSpeciesAny(alternativeChargedSpeciesNames("K+"));   // K+, K[+]
    const Index iCa  = mixture.indexChargedSpeciesAny(alternativeChargedSpeciesNames("Ca++")); // Ca++, Ca+2, Ca[+2]
    const Index iMg  = mixture.indexChargedSpeciesAny(alternativeChargedSpeciesNames("Mg++")); // Mg++, Mg+2, Mg[+2]
    const Index iCl  = mixture.indexChargedSpeciesAny(alternativeChargedSpeciesNames("Cl-"));  // Cl-, Cl[-]

    // The molalities of some ionic species covered by the model
    real mNa(nspecies);
    real mK(nspecies);
    real mCa(nspecies);
    real mMg(nspecies);
    real mCl(nspecies);
    real mSO4(nspecies);

    // The ln activity coefficient of CO2(aq)
    real ln_gCO2(nspecies);

    AqueousActivityModel f = [=](const AqueousMixtureState& state) mutable
    {
        // Extract temperature from the parameters
        const real& T = state.T;

        // The stoichiometric molalities of the ions in the aqueous mixture and their molar derivatives
        const VectorXd& ms = state.ms;

        // Extract the stoichiometric molalities of the specific ions and their molar derivatives
        if(iNa < nions) mNa = ms[iNa];
        if(iK  < nions) mK  = ms[iK];
        if(iCa < nions) mCa = ms[iCa];
        if(iMg < nions) mMg = ms[iMg];
        if(iCl < nions) mCl = ms[iCl];

        // The Pitzer's parameters of the Rumpf et al. (1994) model
        const real B = 0.254 - 76.82/T - 10656.0/(T*T) + 6312.0e+3/(T*T*T);
        const double Gamma = -0.0028;

        // The ln activity coefficient of CO2(aq)
        real ln_gCO2 = 2*B*(mNa + mK + 2*mCa + 2*mMg) + 3*Gamma*(mNa + mK + mCa + mMg)*mCl;

        return ln_gCO2;
    };

    return f;
}

} // namespace Reaktoro
