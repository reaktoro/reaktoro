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

#include "ActivityModelIdealGas.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>

namespace Reaktoro {

using std::log;

auto ActivityModelIdealGas() -> ActivityModelGenerator
{
    ActivityModelGenerator model = [](const SpeciesList& species)
    {
        const auto R = universalGasConstant;

        ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args)
        {
            const auto& [T, P, x] = args;

            const auto Pbar = P * 1.0e-5; // from Pa to bar

            // Set the state of matter of the phase
            props.som = StateOfMatter::Gas;

            props = 0.0;
            props.Vx  =  R*T/P; // identical to entire volume, since V0 = 0 for gases
            props.VxT =  props.Vx/T;
            props.VxP = -props.Vx/P;
            props.ln_a = x.log() + log(Pbar);
        };

        return fn;
    };

    return model;
}

} // namespace Reaktoro
