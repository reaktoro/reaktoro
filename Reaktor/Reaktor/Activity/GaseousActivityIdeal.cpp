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

#include "GaseousActivityIdeal.hpp"

// C++ includes
#include <functional>
using namespace std::placeholders;

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Mixtures/GaseousMixture.hpp>
#include <Reaktor/Common/ConvertUtils.hpp>

namespace Reaktor {
namespace internal {

auto gaseousActivityIdeal(const GaseousActivityParams& params, Index ispecies) -> ThermoScalar
{
    // The pressure (in units of bar)
    const double Pb = convert<Pa,bar>(params.P);

    // The molar fractions of all gaseous species
    const auto& x = params.x;

    // The molar fraction of the given gaseous species
    ThermoScalar xi = x.row(ispecies);

    // The activity of the given gaseous species
    ThermoScalar ai;
    ai.val = xi.val * Pb;
    ai.ddn = xi.ddn * Pb;

    return ai;
}

} /* namespace internal */

auto gaseousActivityIdeal(const std::string& species, const GaseousMixture& mixture) -> GaseousActivity
{
    const Index ispecies = indexSpecies(mixture, species);

    return std::bind(internal::gaseousActivityIdeal, _1, ispecies);
}

} // namespace Reaktor
