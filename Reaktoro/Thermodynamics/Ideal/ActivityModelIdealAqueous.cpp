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

#include "ActivityModelIdealAqueous.hpp"

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {

using std::log;

auto ActivityModelIdealAqueous() -> ActivityModelGenerator
{
    ActivityModelGenerator model = [](const SpeciesList& species)
    {
        const auto iH2O = species.indexWithFormula("H2O");
        const auto MH2O = waterMolarMass;

        ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args)
        {
            const auto x = args.x;
            const auto m = x/(MH2O * args.x[iH2O]); // molalities
            props = 0.0;
            props.ln_a = m.log();
            props.ln_a[iH2O] = log(x[iH2O]);
        };

        return fn;
    };

    return model;
}

} // namespace Reaktoro
