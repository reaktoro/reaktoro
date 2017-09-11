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

#include "MineralMixture.hpp"

namespace Reaktoro {

MineralMixture::MineralMixture()
: GeneralMixture<MineralSpecies>()
{}

MineralMixture::MineralMixture(const std::vector<MineralSpecies>& species)
: GeneralMixture<MineralSpecies>(species)
{}

MineralMixture::MineralMixture(const MineralSpecies& species)
: GeneralMixture<MineralSpecies>({species})
{}

MineralMixture::~MineralMixture()
{}

auto MineralMixture::state(Temperature T, Pressure P, VectorConstRef n) const -> MineralMixtureState
{
    MineralMixtureState res;
    res.T = T;
    res.P = P;
    res.x = molarFractions(n);
    return res;
}

} // namespace Reaktoro
