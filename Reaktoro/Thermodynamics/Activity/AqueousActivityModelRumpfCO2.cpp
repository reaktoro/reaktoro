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

#include "AqueousActivityModelRumpfCO2.hpp"

// Reaktoro includes
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Extensions/Geochemistry/AqueousMixture.hpp>

namespace Reaktoro {

auto aqueousActivityModelRumpfCO2(const AqueousMixture& mixture) -> AqueousActivityModel
{
    // The number of speciesn and charged species
    const auto nspecies = mixture.species().size();
    const auto nions = mixture.charged().size();

    // The local indices of some charged species among all charged species
    const auto iNa  = mixture.charged().indexWithFormula("Na+");
    const auto iK   = mixture.charged().indexWithFormula("K+");
    const auto iCa  = mixture.charged().indexWithFormula("Ca++");
    const auto iMg  = mixture.charged().indexWithFormula("Mg++");
    const auto iCl  = mixture.charged().indexWithFormula("Cl-");

    AqueousActivityModel f = [=](const AqueousMixtureState& state) mutable
    {
        // Extract temperature from the parameters
        const auto T = state.T;

        // The stoichiometric molalities of the ions in the aqueous mixture and their molar derivatives
        const auto& ms = state.ms;

        // Extract the stoichiometric molalities of specific ions
        const real mNa = (iNa < nions) ? ms[iNa] : 0.0;
        const real mK  = (iK  < nions) ? ms[iK]  : 0.0;
        const real mCa = (iCa < nions) ? ms[iCa] : 0.0;
        const real mMg = (iMg < nions) ? ms[iMg] : 0.0;
        const real mCl = (iCl < nions) ? ms[iCl] : 0.0;

        // The Pitzer's parameters of the Rumpf et al. (1994) model
        const auto B = 0.254 - 76.82/T - 10656.0/(T*T) + 6312.0e+3/(T*T*T);
        const auto Gamma = -0.0028;

        // Return the ln activity coefficient of CO2(aq)
        return 2*B*(mNa + mK + 2*mCa + 2*mMg) + 3*Gamma*(mNa + mK + mCa + mMg)*mCl;
    };

    return f;
}

} // namespace Reaktoro
