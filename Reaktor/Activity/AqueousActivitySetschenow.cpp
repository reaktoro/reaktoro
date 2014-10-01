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

#include "AqueousActivitySetschenow.hpp"

// C++ includes
#include <functional>
using namespace std::placeholders;

// Reaktor includes
#include <Reaktor/Common/Index.hpp>

namespace Reaktor {
namespace internal {

auto aqueousActivitySetschenow(const AqueousSolutionState& state, Index ispecies, Index iwater, double b) -> ChemicalScalar
{
    // The effective ionic strength of the aqueous solution
    const auto& I = state.Ie;

    // The molar fractions of the aqueous species and their molar derivatives
    const auto& x = state.x;

    // The molalities of the aqueous species and their molar derivatives
    const auto& m = state.m;

    // The molar fractions of the aqueous species H2O(l) and its molar derivatives
    const double xw_val = x.val().at(iwater);
    const Vector xw_ddn = x.ddn().row(iwater);

    // The molality of the given aqueous species and its molar derivatives
    const double mi_val = m.val().at(ispecies);
    const Vector mi_ddn = m.ddn().row(ispecies);

    // The activity coefficient of the given species and its molar derivatives
    const double gi_val = xw_val * std::pow(10.0, b * I.val());
    const Vector gi_ddn = gi_val * (xw_ddn/xw_val + 2.303*b*I.ddn());

    // The activity of the given species and its molar derivatives
    const double ai_val = mi_val * gi_val;
    const Vector ai_ddn = mi_val * gi_ddn + mi_ddn * gi_val;

    return {ai_val, 0.0, 0.0, ai_ddn};
}

} /* namespace */

auto aqueousActivitySetschenow(const std::string& species, const AqueousSolution& solution, double b) -> AqueousActivity
{
    const Index ispecies = speciesIndex(solution, species);
    const Index iwater   = waterIndex(solution);

    return std::bind(internal::aqueousActivitySetschenow, _1, ispecies, iwater, b);
}

} // namespace Reaktor
