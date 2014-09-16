// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "AqueousActivityIdeal.hpp"

// C++ includes
#include <functional>
using namespace std::placeholders;

// Reaktor includes
#include <Reaktor/Common/Index.hpp>

namespace Reaktor {
namespace internal {

auto aqueousActivityIdeal(const AqueousMixtureState& params, Index ispecies, Index iwater) -> ThermoScalar
{
    // The molar fractions of the aqueous species
    const auto& x = params.x;

    // The molalities of the aqueous species
    const auto& m = params.m;

    // The molar fraction of the aqueous species H2O(l) and its molar derivatives
    ThermoScalar xw = x.row(iwater);

    // The molality of the given aqueous species and its molar derivatives
    ThermoScalar mi = m.row(ispecies);

    // The activity of the given aqueous species and its molar derivatives
    ThermoScalar ai;
    ai.val = mi.val * xw.val;
    ai.ddn = mi.val * xw.ddn + mi.ddn * xw.val;

    return ai;
}

auto aqueousActivityIdealWater(const AqueousMixtureState& params, Index iwater) -> ThermoScalar
{
    // The molar fractions of the aqueous species
    const auto& x = params.x;

    // The molar fraction of the aqueous species H2O(l) and its molar derivatives
    ThermoScalar xw = x.row(iwater);

    return xw;
}

} /* namespace internal */

auto aqueousActivityIdeal(const std::string& species, const AqueousMixture& mixture) -> AqueousActivity
{
    const Index ispecies = indexSpecies(mixture, species);
    const Index iwater = indexWater(mixture);

    if(ispecies == iwater) return std::bind(internal::aqueousActivityIdealWater, _1, iwater);
    else return std::bind(internal::aqueousActivityIdeal, _1, ispecies, iwater);
}

} // namespace Reaktor
