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

#include "ActivityModelIdealAqueous.hpp"

// Reaktoro includes
#include <Reaktoro/Water/WaterConstants.hpp>

namespace Reaktoro {

using std::log;

auto ActivityModelIdealAqueous() -> ActivityModelGenerator
{
    ActivityModelGenerator model = [](const SpeciesList& species)
    {
        const auto iw = species.indexWithFormula("H2O");
        const auto Mw = waterMolarMass;

        ActivityModel fn = [=](ActivityPropsRef props, ActivityModelArgs args)
        {
            const auto x = args.x;
            const auto xw = x[iw];
            const auto m = x/(Mw * xw); // molalities

            // Set the state of matter of the phase
            props.som = StateOfMatter::Liquid;

            props = 0.0;
            props.ln_a = m.log();
            props.ln_a[iw] = -(1 - xw)/xw; // consistent to Gibbs-Duhem conditions
        };

        return fn;
    };

    return model;
}

} // namespace Reaktoro
