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
namespace {

auto computeAqueousActivityIdeal(const AqueousMixtureState& state, Index ispecies, Index iwater) -> ChemicalScalar
{
    // The molar fractions of the aqueous species
    const auto& x = state.x;

    // The molalities of the aqueous species
    const auto& m = state.m;

    // The molar fraction of the aqueous species H2O(l) and its molar derivatives
    ChemicalScalar xw = x.row(iwater);

    // The molality of the given aqueous species and its molar derivatives
    ChemicalScalar mi = m.row(ispecies);

    // The activity of the given aqueous species and its molar derivatives
    ChemicalScalar ai;
    ai.val = mi.val * xw.val;
    ai.ddn = mi.val * xw.ddn + mi.ddn * xw.val;

    return ai;
}

auto computeAqueousActivityIdealWater(const AqueousMixtureState& state, Index iwater) -> ChemicalScalar
{
    // The molar fractions of the aqueous species
    const auto& x = state.x;

    // The molar fraction of the aqueous species H2O(l) and its molar derivatives
    ChemicalScalar xw = x.row(iwater);

    return xw;
}

} // namespace

auto aqueousActivityIdeal(const std::string& species, const AqueousMixture& mixture) -> AqueousActivityFunction
{
    const Index ispecies = mixture.indexSpecies(species);
    const Index iwater = mixture.indexSpecies("H2O(l)");

    if(ispecies == iwater) return std::bind(computeAqueousActivityIdealWater, _1, iwater);
    else return std::bind(computeAqueousActivityIdeal, _1, ispecies, iwater);
}

} // namespace Reaktor
