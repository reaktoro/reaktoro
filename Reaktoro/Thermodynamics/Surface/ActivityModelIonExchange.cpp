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
#include <Reaktoro/Extensions/Phreeqc/PhreeqcLegacy.hpp>

namespace Reaktoro {

using std::sqrt;
using std::log;

namespace detail {

// Return the number of exchanger's equivalents (the charge of cations) in the ion exchange species.
auto exchangerEquivalentsNumber(const Species& species) -> real
{
    // Run through the elements of the current species and return the coefficient of the exchanger
    for(auto [element, coeff] : species.elements())
        if(!Elements::withSymbol(element.symbol()))
            return coeff;

    // If all the elements are part of the periodic table then the exchanger is missing
    errorif(true, "Could not get information about the exchanger equivalents number. "
        "Ensure the ion exchange phase contains correct species");
}

/// Return the IonExchangeActivityModel object based on the Gaines--Thomas model.
auto activityModelIonExchangeGainesThomas(const SpeciesList& species) -> ActivityModel
{
    // The number of species in the ion exchange phase only
    const auto num_species = species.size();

    // The numbers of exchanger's equivalents for exchange species
    ArrayXd ze = ArrayXr::Zero(num_species);

    // Initialize exchanger's equivalents by parsing the elements of the ion exchange species
    for(auto i = 0; i < num_species; ++i)
        ze[i] = detail::exchangerEquivalentsNumber(species[i]);

    // Define the activity model function of the ion exchange phase
    ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x] = args;

        // Auxiliary references
        auto& ln_g = props.ln_g;
        auto& ln_a = props.ln_a;

        // Export the exchanger equivalents of ion exchange composition via `extra` data member
        props.extra["ExchangerEquivalents"] = ze;

        // Calculate the ln of equivalence fractions
        const auto ln_beta = (x*ze/(x*ze).sum()).log();

        // Calculate the ln of activity coefficients
        ln_g = ArrayXr::Zero(num_species);

        // Calculate Davies and Debye--Huckel parameters only if the AqueousPhase has been already evaluated
        if (props.extra["AqueousMixtureState"].has_value())
        {
            // Export aqueous mixture state via `extra` data member
            const auto& state = std::any_cast<AqueousMixtureState>(props.extra["AqueousMixtureState"]);

            // Auxiliary constant references properties
            const auto& I = state.Is;            // the stoichiometric ionic strength
            const auto& rho = state.rho/1000;    // the density of water (in g/cm3)
            const auto& epsilon = state.epsilon; // the dielectric constant of water

            // Auxiliary variables
            const auto sqrtI = sqrt(I);
            const auto sqrt_rho = sqrt(rho);
            const auto T_epsilon = T * epsilon;
            const auto sqrt_T_epsilon = sqrt(T_epsilon);
            const auto A = 1.824829238e+6 * sqrt_rho/(T_epsilon*sqrt_T_epsilon);
            const auto B = 50.29158649 * sqrt_rho/sqrt_T_epsilon;
            const auto ln10 = log(10);

            // Loop over all species in the composition
            for(Index i = 0; i < num_species; ++i)
            {
                //            const auto phreeqc_species = std::any_cast<PhreeqcSpecies*>(species[i].attachedData());
                //            std::cout << phreeqc_species->name << "" << phreeqc_species->dw << std::endl;
                //            std::cout << phreeqc_species->dha << std::endl;
                //            std::cout << phreeqc_species->dhb << std::endl;
                //            std::cout << phreeqc_species->a_f << std::endl;
                //            getchar();
                // Calculate activity coefficients according to the Debye--Huckel model
                // TODO: obtained from each species if it have parameter -gamma provided
                //            const auto a = phreeqc_species->dha;
                //            const auto b = phreeqc_species->dhb;
                const auto a = 1.0;
                const auto b = 1.0;

                // Calculate the ln activity coefficient of the exchange species
                ln_g[i] = ln10*(-A*ze[i]*ze[i]*sqrtI/(1.0 + a*B*sqrtI) + b*I);

                // ---------------------------------------------------------------------------//
                // Calculate activity coefficients according top the Davies model

                // Calculate the ln activity coefficient of the echange species
                // Debye-Huckel parameter
                const auto Agamma = 0.5095;
                ln_g[i] = ln10*(-Agamma*ze[i]*ze[i]*sqrtI/(1 + sqrtI) - 0.3 * I);
            }
        }

        // Calculate the ln of activities
        ln_a = ln_g + ln_beta;
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
