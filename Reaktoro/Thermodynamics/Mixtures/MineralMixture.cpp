// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

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

auto MineralMixture::state(Temperature T, Pressure P, VectorXrConstRef n) const -> MineralMixtureState
{
    MineralMixtureState res;
    res.T = T;
    res.P = P;
    res.x = moleFractions(n);
    return res;
}

} // namespace Reaktoro
