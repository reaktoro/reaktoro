// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

#include "GaseousMixture.hpp"

namespace Reaktoro {

GaseousMixture::GaseousMixture()
: GeneralMixture<GaseousSpecies>()
{}

GaseousMixture::GaseousMixture(const std::vector<GaseousSpecies>& species)
: GeneralMixture<GaseousSpecies>(species)
{}

GaseousMixture::~GaseousMixture()
{}

auto GaseousMixture::state(Temperature T, Pressure P, VectorConstRef n) const -> GaseousMixtureState
{
    GaseousMixtureState res;
    res.T = T;
    res.P = P;
    res.x = molarFractions(n);
    return res;
}

} // namespace Reaktoro
