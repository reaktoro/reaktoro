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

#include "FluidMixture.hpp"
#include <Reaktoro/Util/PhaseIdentificationMethods.hpp>

namespace Reaktoro {

FluidMixture::FluidMixture()
    : GeneralMixture<FluidSpecies>()
{}

FluidMixture::FluidMixture(const std::vector<FluidSpecies>& species)
    : GeneralMixture<FluidSpecies>(species)
{}

FluidMixture::FluidMixture(const std::vector<FluidSpecies>& species, PhaseID::PhaseIdentificationMethods fluidMixturePhaseIdentificationMethod)
    : GeneralMixture<FluidSpecies>(species), _fluidMixturePhaseIdentificationMethod(fluidMixturePhaseIdentificationMethod)
{}

FluidMixture::~FluidMixture()
{}

auto FluidMixture::fluidMixturePhaseIdentificationMethod() const -> PhaseID::PhaseIdentificationMethods
{
    return _fluidMixturePhaseIdentificationMethod;
}

auto FluidMixture::setFluidMixturePhaseIdentificationMethod(PhaseID::PhaseIdentificationMethods fluidMixturePhaseIdentificationMethod) -> void
{
    _fluidMixturePhaseIdentificationMethod = fluidMixturePhaseIdentificationMethod;
}

auto FluidMixture::removeInapproprieatePhase() const -> bool
{
    return _RemoveInapproprieatePhase;
}

auto FluidMixture::setRemoveInapproprieatePhaseAsTrue() -> void
{
    _RemoveInapproprieatePhase = true;
}

auto FluidMixture::setRemoveInapproprieatePhaseAsFalse() -> void
{
    _RemoveInapproprieatePhase = false;
}

auto FluidMixture::state(Temperature T, Pressure P, VectorConstRef n) const -> FluidMixtureState
{
    FluidMixtureState res;
    res.T = T;
    res.P = P;
    res.x = moleFractions(n);
    return res;
}

} // namespace Reaktoro
