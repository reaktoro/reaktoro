// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

#include "ActivityModelIonExchange.hpp"

// Reaktoro includes
#include <Reaktoro/Singletons/Elements.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousProps.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousMixture.hpp>
#include <Reaktoro/Thermodynamics/Surface/IonExchangeSurface.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcLegacy.hpp>

namespace Reaktoro {

using std::sqrt;
using std::log;

namespace detail {

/// Return the IonExchangeActivityModel object based on the Gaines--Thomas model.
auto activityModelIonExchangeGainesThomas(const SpeciesList& species) -> ActivityModel
{
    // Create the ion exchange surface
    IonExchangeSurface surface(species);

    // The number of all species and ion exchange species in the current exchange phase only
    const auto num_species = species.size();

    // The numbers of exchanger's equivalents for exchange species
    ArrayXd ze = surface.ze();

    // Define the activity model function of the ion exchange phase
    ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x] = args;

        // Auxiliary references
        auto& ln_g = props.ln_g;
        auto& ln_a = props.ln_a;

        // Calculate ln of activities of ion exchange species as the ln of equivalence fractions
        ln_a = (x*ze/(x*ze).sum()).log();

        // Initialized the ln of activity coefficients of the ion exchange species on the surface
        ln_g = ArrayXr::Zero(num_species);

        // Calculate Davies and Debye--Huckel parameters only if the AqueousPhase has been already evaluated
        if (props.extra["AqueousMixtureState"].has_value())
        {
            // Export aqueous mixture state via `extra` data member
            const auto& aqstate = std::any_cast<AqueousMixtureState>(props.extra["AqueousMixtureState"]);

            // Auxiliary constant references properties
            const auto& I = aqstate.Is;            // the stoichiometric ionic strength
            const auto& rho = aqstate.rho/1000;    // the density of water (in g/cm3)
            const auto& epsilon = aqstate.epsilon; // the dielectric constant of water

            // Auxiliary variables
            const auto sqrtI = sqrt(I);
            const auto sqrt_rho = sqrt(rho);
            const auto T_epsilon = T * epsilon;
            const auto sqrt_T_epsilon = sqrt(T_epsilon);
            const auto A = 1.824829238e+6 * sqrt_rho/(T_epsilon*sqrt_T_epsilon);
            const auto B = 50.29158649 * sqrt_rho/sqrt_T_epsilon;
            const auto ln10 = log(10);
            const auto Agamma = 0.5095; // the Debye-Huckel parameter

            // Loop over all ion exchange species in composition
            for(auto i = 0; i < num_species; i++)
            {
                // Fetch phreeqc species from the `attachedData` field of the species
                const auto phreeqc_species = std::any_cast<const PhreeqcSpecies*>(species[i].attachedData());

                // Fetch Debye--Huckel activity model parameters a and b
                const auto a = phreeqc_species->dha;
                const auto b = phreeqc_species->dhb;

                if(b == 99.9) // the Phreeqc hack used while reading the dat-file, which indicate to use Davies activity model
                {
                    // Calculate the ln activity coefficient of the exchange species using the Davies activity model
                    ln_g[i] = ln10*(-Agamma*ze[i]*ze[i]*sqrtI/(1 + sqrtI) - 0.3*I);
                }
                else // otherwise, the Debye--Huckel activity model is preferred
                {
                    // Calculate the ln activity coefficient of the exchange species using the Debye--Huckel model
                    ln_g[i] = ln10*(-A*ze[i]*ze[i]*sqrtI/(1.0 + a*B*sqrtI) + b*I);
                }
            }
        }
        // Add the correction introduced by the activity coefficients
        ln_a += ln_g;

    };
    return fn;
}

} // namespace detail

auto ActivityModelIonExchange() -> ActivityModelGenerator
{
    return ActivityModelIonExchangeGainesThomas();
}

auto ActivityModelIonExchangeGainesThomas() -> ActivityModelGenerator
{
    return [=](const SpeciesList& species)
    {
        return detail::activityModelIonExchangeGainesThomas(species);
    };
}

} // namespace Reaktoro
