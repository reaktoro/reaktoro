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

#pragma once

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Species/GaseousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/GeneralMixture.hpp>

namespace Reaktoro {

/// A type used to describe the state of a gaseous mixture
struct GaseousMixtureState : public MixtureState
{};

/// Provides a computational representation of a gaseous mixture.
/// The GaseousMixture class is defined as a collection of GaseousSpecies objects,
/// representing, therefore, a mixture of gaseous species. Its main purpose is to
/// provide the necessary operations in the calculation of activities of gaseous
/// species.
/// @see GaseousSpecies
/// @ingroup Mixtures
class GaseousMixture : public GeneralMixture<GaseousSpecies>
{
public:
    /// Construct a default GaseousMixture instance.
    GaseousMixture();

    /// Construct a GaseousMixture instance with given species.
    /// @param species The species that compose the gaseous mixture
    explicit GaseousMixture(const std::vector<GaseousSpecies>& species);

    /// Destroy the GaseousMixture instance.
    virtual ~GaseousMixture();

    /// Calculate the state of the gaseous mixture.
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    /// @param n The molar amounts of the species in the mixture (in units of mol)
    auto state(Temperature T, Pressure P, const Vector& n) const -> GaseousMixtureState;
};

} // namespace Reaktoro
