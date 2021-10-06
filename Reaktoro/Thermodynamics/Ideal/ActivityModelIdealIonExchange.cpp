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
#include <Reaktoro/Thermodynamics/Surface/IonExchangeSurface.hpp>

namespace Reaktoro {

auto ActivityModelIdealIonExchange() -> ActivityModelGenerator
{
    ActivityModelGenerator model = [](const SpeciesList& species)
    {
        // Create the ion exchange surface
        IonExchangeSurface surface(species);

        // The numbers of exchanger's equivalents for exchange species
        ArrayXd ze = surface.ze();

        ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args)
        {
            // Fetch species fractions for the activity model evaluation
            const auto x = args.x;

            // Calculate ln of activities of ion exchange species as the ln of equivalence fractions
            props.ln_a = (x*ze/(x*ze).sum()).log();
        };

        return fn;
    };

    return model;
}

} // namespace Reaktoro
