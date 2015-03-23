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

#include "MineralActivityIdeal.hpp"

// C++ includes
#include <functional>
using namespace std::placeholders;

// Reaktor includes
#include <Reaktor/Common/Index.hpp>

namespace Reaktor {
namespace {

auto computeMineralActivityIdeal(const MineralMixtureState& state, Index ispecies) -> ChemicalScalar
{
    return state.x.row(ispecies);
}

} // namespace

auto mineralActivityIdeal(const std::string& species, const MineralMixture& mixture) -> MineralActivityFunction
{
    const Index ispecies = mixture.indexSpecies(species);

    return std::bind(computeMineralActivityIdeal, _1, ispecies);
}

} // namespace Reaktor
