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

#include "ActivityModelPhreeqcIonicStrengthPressureCorrection.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcLegacy.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcThermo.hpp>
#include <Reaktoro/Models/ActivityModels/Support/AqueousMixture.hpp>

namespace Reaktoro {

auto ActivityModelPhreeqcIonicStrengthPressureCorrection() -> ActivityModelGenerator
{
    ActivityModelGenerator model = [=](const SpeciesList& specieslist)
    {
        // Collect the PhreeqcSpecies pointers attached to each Species object
        Vec<const PhreeqcSpecies*> phspecies;
        for(auto species : specieslist)
        {
            const auto has_phreeqc_species = std::any_cast<const PhreeqcSpecies*>(&species.attachedData());
            if(has_phreeqc_species)
                phspecies.push_back(std::any_cast<const PhreeqcSpecies*>(species.attachedData()));
            else phspecies.push_back(nullptr); // if PhreeqcSpecies pointer not attached, consider null pointer
        }

        const auto Pref = 1.0e5; // reference pressure at 1 bar (in Pa)

        ActivityModel fn = [=](ActivityPropsRef props, ActivityModelArgs args)
        {
            // The arguments for the activity model evaluation
            const auto& [T, P, x] = args;

            // Check AqueousMixtureState is available in props.extra
            auto stateit = props.extra.find("AqueousMixtureState");

            errorif(stateit == props.extra.end(),
                "ActivityModelPhreeqcIonicStrengthPressureCorrection expects that another aqueous activity model has been chained first (e.g., Davies, Debye-Huckel, HKF, PitzerHMW, etc.) ");

            // The aqueous mixture state exported by a base aqueous activity model.
            const auto& state = *std::any_cast<SharedPtr<AqueousMixtureState> const&>(stateit->second);

            const auto mu = state.Ie;
            const auto RT = universalGasConstant * T;

            const auto wprops = PhreeqcUtils::waterPropsMemoized(T, P);

            for(auto [i, species] : enumerate(specieslist))
            {
                if(phspecies[i] == nullptr)
                    continue;

                const auto Vcorr = standardVolumeIonicStrengthCorrection(phspecies[i], T, P, mu, wprops) * cubicCentimeterToCubicMeter; // in m3/mol

                props.ln_a[i] += Vcorr * (P - Pref)/RT;
            }
        };

        return fn;
    };

    return model;
}

} // namespace Reaktoro
