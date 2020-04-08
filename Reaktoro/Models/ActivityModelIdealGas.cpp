// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

auto ActivityModelIdealGas::create(const SpeciesList& species) -> ActivityModelFn
{
    ActivityModelFn fn = [](ActivityProps props, real T, real P, ArrayXrConstRef x)
    {
        using std::log;
        const auto R = universalGasConstant;
        const auto Pbar = P * 1.0e-5; // from Pa to bar
        props.Vex  = 0.0;
        props.VexT = 0.0;
        props.VexP = 0.0;
        props.Gex  = 0.0;
        props.Hex  = 0.0;
        props.Cpex = 0.0;
        props.ln_g = 0.0;
        props.ln_a = x.log() + log(Pbar);
    };

    return fn;
}

} // namespace Reaktoro
