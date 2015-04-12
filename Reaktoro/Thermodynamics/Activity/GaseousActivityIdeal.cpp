// Reaktoro is a C++ library for computational reaction modelling.
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

#include "GaseousActivityIdeal.hpp"

// C++ includes
#include <functional>
using namespace std::placeholders;

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/ConvertUtils.hpp>

namespace Reaktoro {
namespace {

auto computeGaseousActivityIdeal(const GaseousMixtureState& state, Index ispecies) -> ChemicalScalar
{
    // The pressure (in units of bar)
    const double Pb = convert<Pa,bar>(state.P);

    // The molar fractions of all gaseous species
    const auto& x = state.x;

    // The molar fraction of the given gaseous species
    ChemicalScalar xi = x.row(ispecies);

    // The activity of the given gaseous species
    ChemicalScalar ai;
    ai.val = xi.val * Pb;
    ai.ddn = xi.ddn * Pb;

    return ai;
}

} // namespace

auto gaseousActivityIdeal(const std::string& species, const GaseousMixture& mixture) -> GaseousActivityFunction
{
    const Index ispecies = mixture.indexSpecies(species);

    return std::bind(computeGaseousActivityIdeal, _1, ispecies);
}

} // namespace Reaktoro
