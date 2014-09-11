/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "AqueousActivitySetschenow.hpp"

// C++ includes
#include <functional>
using namespace std::placeholders;

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Mixtures/AqueousMixture.hpp>

namespace Reaktor {
namespace internal {

auto aqueousActivitySetschenow(const AqueousActivityParams& params, Index ispecies, Index iwater, double b) -> ThermoScalar
{
    // The effective ionic strength of the aqueous mixture
    const auto& I = params.Ie;

    // The molar fractions of the aqueous species and their molar derivatives
    const auto& x = params.x;

    // The molalities of the aqueous species and their molar derivatives
    const auto& m = params.m;

    // The molar fractions of the aqueous species H2O(l) and its molar derivatives
    const ThermoScalar xw = x.row(iwater);

    // The molality of the given aqueous species and its molar derivatives
    const ThermoScalar mi = m.row(ispecies);

    // The activity coefficient of the given species and its molar derivatives
    ThermoScalar gi;
    gi.val = xw.val * std::pow(10.0, b * I.val);
    gi.ddn = gi.val * (xw.ddn/xw.val + 2.303*b*I.ddn);

    // The activity of the given species and its molar derivatives
    ThermoScalar ai;
    ai.val = mi.val * gi.val;
    ai.ddn = mi.val * gi.ddn + mi.ddn * gi.val;

    return ai;
}

} /* namespace */

auto aqueousActivitySetschenow(const std::string& species, const AqueousMixture& mixture, double b) -> AqueousActivity
{
    const Index ispecies = indexSpecies(mixture, species);
    const Index iwater   = mixture.indexWater();

    return std::bind(internal::aqueousActivitySetschenow, _1, ispecies, iwater, b);
}

} // namespace Reaktor
