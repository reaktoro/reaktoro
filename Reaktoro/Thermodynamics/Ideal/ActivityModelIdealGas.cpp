// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

namespace Reaktoro {

using std::log;

auto ActivityModelIdealGas::operator()(const SpeciesList& species) -> ActivityPropsFn
{
    ActivityPropsFn fn = [](ActivityProps props, ActivityArgs args) mutable
    {
        const auto Pbar = args.P * 1.0e-5; // from Pa to bar
        props = 0.0;
        props.ln_a = args.x.log() + log(Pbar);
    };

    return fn;
}

} // namespace Reaktoro
