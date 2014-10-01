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

#include "GaseousActivityIdeal.hpp"

// C++ includes
#include <functional>
using namespace std::placeholders;

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/ConvertUtils.hpp>

namespace Reaktor {
namespace internal {

auto gaseousActivityIdeal(const GaseousSolutionState& state, Index ispecies) -> ChemicalScalar
{
    // The pressure (in units of bar)
    const double Pb = convert<Pa,bar>(state.P);

    // The molar fractions of all gaseous species
    const auto& x = state.x;

    // The molar fraction of the given gaseous species
    const double xi_val = x.val().at(ispecies);
    const Vector xi_ddn = x.ddn().row(ispecies);

    // The activity of the given gaseous species
    const double ai_val = xi_val * Pb;
    const Vector ai_ddn = xi_ddn * Pb;

    return {ai_val, 0.0, 0.0, ai_ddn};
}

} /* namespace internal */

auto gaseousActivityIdeal(const std::string& species, const GaseousSolution& solution) -> GaseousActivity
{
    const Index ispecies = speciesIndex(solution, species);

    return std::bind(internal::gaseousActivityIdeal, _1, ispecies);
}

} // namespace Reaktor
