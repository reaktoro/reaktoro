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

namespace Reaktoro {

using std::sqrt;

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
/// Return the InoExchangeActivityModel object based on the Gaines--Thomas model.
auto activityModelIonExchangeGainesThomas(const SpeciesList& species) -> ActivityModel
{
    // The number of species
    auto num_species = species.size();

    // The numbers of exchanger's equivalents for exchange species
    ArrayXd ze = ArrayXr::Zero(num_species);
    // Initialize exchanger's equivalents by parsing the elements of the ion exchange species
    for(Index i = 0; i < num_species; ++i)
        ze[i] = detail::exchangerEquivalentsNumber(species[i]);

    // Define the activity model function of the ion exchange phase
    ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x] = args;

        // Auxiliary references
        auto& ln_g = props.ln_g;
        auto& ln_a = props.ln_a;

        // Export the aqueous mixture and its state via the `extra` data member
        props.extra = { ze };

        // Calculate the ln of equivalence fractions
        ArrayXr ln_beta = (x * ze / (x * ze).sum()).log();

        // Calculate the ln of activity coefficients
        ln_g = ArrayXr::Zero(num_species);

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
