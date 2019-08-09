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
#include <Reaktoro/Util/PhaseIdentificationMethods.hpp>

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
    explicit FluidMixture(const std::vector<FluidSpecies>& species, PhaseID::PhaseIdentificationMethods fluidMixturePhaseIdentificationMethod);

    /// Destroy the FluidMixture instance.
    virtual ~FluidMixture();

    /// Return phase identification method used by the mixture
    auto fluidMixturePhaseIdentificationMethod() const->PhaseID::PhaseIdentificationMethods;

    /// Return if the mixture should remove inappropriate phase
    auto removeInapproprieatePhase() const -> bool;

    /// Set the _RemoveInapproprieatePhase as true
    auto setRemoveInapproprieatePhaseAsTrue() -> void;

    /// Set the _RemoveInapproprieatePhase as false
    auto setRemoveInapproprieatePhaseAsFalse() -> void;

    /// Set the phase identification method used by the mixture
    auto setFluidMixturePhaseIdentificationMethod(PhaseID::PhaseIdentificationMethods  fluidMixturePhaseIdentificationMethod) -> void;


    /// Calculate the state of the fluid (gaseous or liquid) mixture.
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    /// @param n The molar amounts of the species in the mixture (in units of mol)
    auto state(Temperature T, Pressure P, VectorConstRef n) const->FluidMixtureState;

private:
        
    /// The phase identification method for the fluid mixture
    PhaseID::PhaseIdentificationMethods _fluidMixturePhaseIdentificationMethod = PhaseID::PhaseIdentificationMethods::GibbsEnergyAndEquationOfStateMethod;

    /// flag that indicate if the CubicEOS can remove any inappropriate phase
    bool _RemoveInapproprieatePhase = false;

};

} // namespace Reaktoro
