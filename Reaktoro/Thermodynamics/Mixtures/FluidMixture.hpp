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

#pragma once

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Species/FluidSpecies.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/GeneralMixture.hpp>
#include <Reaktoro/Util/PhaseIdentification.hpp>

namespace Reaktoro {

/// A type used to describe the state of a fluid (gaseous or liquid) mixture
struct FluidMixtureState : public MixtureState
{};

/// Provides a computational representation of a fluid (gaseous or liquid) mixture.
/// The FluidMixture class is defined as a collection of FluidSpecies objects,
/// representing, therefore, a mixture of fluid (gaseous or liquid) species. Its main purpose is to
/// provide the necessary operations in the calculation of activities of fluid (gaseous or liquid)
/// species.
/// @see FluidSpecies
/// @ingroup Mixtures
class FluidMixture : public GeneralMixture<FluidSpecies>
{
public:
    /// Construct a default FluidMixture instance.
    FluidMixture();

    /// Construct a FluidMixture instance with given species.
    /// @param species The species that compose the fluid (gaseous or liquid) mixture
    explicit FluidMixture(const std::vector<FluidSpecies>& species);

    /// Construct a FluidMixture instance with given species, phase type mixture and phase identification method
    /// @param species The species that compose the fluid (gaseous or liquid) mixture
    /// @param fluidMixturePhaseIdentificationMethod The phase identification method used by that composed mixture
    explicit FluidMixture(const std::vector<FluidSpecies>& species, PhaseIdentificationMethod fluidMixturePhaseIdentificationMethod);

    /// Destroy the FluidMixture instance.
    virtual ~FluidMixture();

    /// Calculate the state of the fluid (gaseous or liquid) mixture.
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    /// @param n The molar amounts of the species in the mixture (in units of mol)
    auto state(Temperature T, Pressure P, VectorConstRef n) const->FluidMixtureState;
};

} // namespace Reaktoro
