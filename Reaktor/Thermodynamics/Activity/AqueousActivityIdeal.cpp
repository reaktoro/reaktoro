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

auto computeAqueousActivityIdeal(const AqueousSolutionState& state, Index ispecies, Index iwater) -> ChemicalScalar
{
    // The molar fractions of the aqueous species
    const auto& x = state.x;

    // The molalities of the aqueous species
    const auto& m = state.m;

    // The molar fraction of the aqueous species H2O(l) and its molar derivatives
    const double xw_val = x.val()[iwater];
    const Vector xw_ddn = x.ddn().row(iwater);

    // The molality of the given aqueous species and its molar derivatives
    const double mi_val = m.val()[ispecies];
    const Vector mi_ddn = m.ddn().row(ispecies);

    // The activity of the given aqueous species and its molar derivatives
    const double ai_val = mi_val * xw_val;
    const Vector ai_ddn = mi_val * xw_ddn + mi_ddn * xw_val;

    return {ai_val, 0.0, 0.0, ai_ddn};
}

auto computeAqueousActivityIdealWater(const AqueousSolutionState& state, Index iwater) -> ChemicalScalar
{
    // The molar fractions of the aqueous species
    const auto& x = state.x;

    // The molar fraction of the aqueous species H2O(l) and its molar derivatives
    const double xw_val = x.val()[iwater];
    const Vector xw_ddn = x.ddn().row(iwater);

    return {xw_val, 0.0, 0.0, xw_ddn};
}

} // namespace

auto aqueousActivityIdeal(const std::string& species, const AqueousSolution& solution) -> AqueousActivity
{
    const Index ispecies = speciesIndex(solution, species);
    const Index iwater = waterIndex(solution);

    if(ispecies == iwater) return std::bind(computeAqueousActivityIdealWater, _1, iwater);
    else return std::bind(computeAqueousActivityIdeal, _1, ispecies, iwater);
}

} // namespace Reaktor
