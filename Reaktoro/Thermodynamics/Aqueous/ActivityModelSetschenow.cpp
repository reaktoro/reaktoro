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

#include "ActivityModelSetschenow.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousMixture.hpp>

namespace Reaktoro {

using std::log;

auto ActivityModelSetschenow(String neutral, real b) -> ActivityModelGenerator
{
    ActivityModelGenerator model = [=](const SpeciesList& species)
    {
        // The index of the neutral aqueous species in the aqueous phase.
        const auto ineutral = species.indexWithFormula(neutral);

        ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args)
        {
            // The aqueous mixture and its state exported by a base aqueous activity model.
            const auto& state = std::any_cast<AqueousMixtureState>(props.extra["AqueousMixtureState"]);

            const auto& I = state.Is;
            props.ln_g[ineutral] = ln10 * b * I;
            props.ln_a[ineutral] = props.ln_g[ineutral] + log(state.m[ineutral]);
        };

        return fn;
    };

    return model;
}

} // namespace Reaktoro
