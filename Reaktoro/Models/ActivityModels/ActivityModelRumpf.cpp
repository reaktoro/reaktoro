// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

#include "ActivityModelRumpf.hpp"

// Reaktoro includes
#include <Reaktoro/Models/ActivityModels/Support/AqueousMixture.hpp>

namespace Reaktoro {

using std::log;

auto ActivityModelRumpf(String gas) -> ActivityModelGenerator
{
    ActivityModelGenerator model = [=](const SpeciesList& species)
    {
        // The index of the dissolved gas in the aqueous phase.
        const auto igas = species.indexWithFormula(gas);

        // the local indices of some charged species among all charged species
        const auto aqmix = AqueousMixture(species);
        const auto iNa  = aqmix.charged().findWithFormula("Na+");
        const auto iK   = aqmix.charged().findWithFormula("K+");
        const auto iCa  = aqmix.charged().findWithFormula("Ca++");
        const auto iMg  = aqmix.charged().findWithFormula("Mg++");
        const auto iCl  = aqmix.charged().findWithFormula("Cl-");

        ActivityModel fn = [=](ActivityPropsRef props, ActivityModelArgs args)
        {
            // Check AqueousMixture and AqueousMixtureState are available in props.extra
            auto mixtureit = props.extra.find("AqueousMixture");
            auto stateit = props.extra.find("AqueousMixtureState");

            errorif(stateit == props.extra.end(),
                "ActivityModelRumpf expects that another aqueous activity model has been chained first (e.g., Davies, Debye-Huckel, HKF, PitzerHMW, etc.) ");

            // The aqueous mixture and its state exported by a base aqueous activity model.
            const auto& mixture = *std::any_cast<SharedPtr<AqueousMixture> const&>(mixtureit->second);
            const auto& state = *std::any_cast<SharedPtr<AqueousMixtureState> const&>(stateit->second);

            // The number of charged species
            const auto nions = mixture.charged().size();

            // Extract temperature from the parameters
            const auto T = state.T;

            // The stoichiometric molalities of the ions in the aqueous mixture and their molar derivatives
            const auto& ms = state.ms;

            // Extract the stoichiometric molalities of specific ions
            const auto mNa = (iNa < nions) ? ms[iNa] : real(0.0);
            const auto mK  = (iK  < nions) ? ms[iK]  : real(0.0);
            const auto mCa = (iCa < nions) ? ms[iCa] : real(0.0);
            const auto mMg = (iMg < nions) ? ms[iMg] : real(0.0);
            const auto mCl = (iCl < nions) ? ms[iCl] : real(0.0);

            // The Pitzer's parameters of the Rumpf et al. (1994) model
            const auto B = 0.254 - 76.82/T - 10656.0/(T*T) + 6312.0e+3/(T*T*T);
            const auto Gamma = -0.0028;

            props.ln_g[igas] = 2*B*(mNa + mK + 2*mCa + 2*mMg) + 3*Gamma*(mNa + mK + mCa + mMg)*mCl;
            props.ln_a[igas] = props.ln_g[igas] + log(state.m[igas]);
        };

        return fn;
    };

    return model;
}

} // namespace Reaktoro
