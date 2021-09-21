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

#include "ActivityModelIdealIonExchange.hpp"

// Reaktoro includes
#include <Reaktoro/Singletons/Elements.hpp>

namespace Reaktoro {

auto ActivityModelIdealIonExchange() -> ActivityModelGenerator
{
    ActivityModelGenerator model = [](const SpeciesList& species)
    {
        // Initialize exchanger's equivalents by parsing the elements of the ion exchange species
        const auto num_species = species.size();
        ArrayXd ze = ArrayXr::Zero(num_species);
        // TODO: replace it with the unified function similar to `detail::exchangerEquivalentsNumber` in `ActivityModelIonExchange.cpp` file
        for(auto i = 0; i < num_species; ++i)
            for(auto [element, coeff] : species[i].elements())
                if(!Elements::withSymbol(element.symbol()))
                    ze[i] = coeff;

        ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args)
        {
            // Fetch species fractions for the activity model evaluation
            const auto x = args.x;

            // Export the exchanger equivalents of ion exchange composition via `extra` data member
            props.extra["ExchangerEquivalents"] = ze;

            // Calculate the ln of activities as lon of equivalence fractions
            props.ln_a = (x*ze/(x*ze).sum()).log();
        };

        return fn;
    };

    return model;
}

} // namespace Reaktoro
