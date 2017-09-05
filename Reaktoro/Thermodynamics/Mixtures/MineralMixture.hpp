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
#include <Reaktoro/Thermodynamics/Species/MineralSpecies.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/GeneralMixture.hpp>

namespace Reaktoro {

/// A type used to describe the state of a mineral mixture
struct MineralMixtureState : public MixtureState
{};

/// Provide a computational representation of a mineral mixture.
/// The MineralMixture class is defined as a collection of MineralSpecies objects,
/// representing, therefore, a mixture of mineral species. Its main purpose is to
/// provide the necessary operations in the calculation of activities of mineral
/// species.
/// @see MineralSpecies
/// @ingroup Mixtures
class MineralMixture : public GeneralMixture<MineralSpecies>
{
public:
    /// Construct a default MineralMixture instance.
    MineralMixture();

    /// Construct a MineralMixture instance with given species.
    /// @param species The species that compose the mineral mixture
    explicit MineralMixture(const std::vector<MineralSpecies>& species);

    /// Construct a MineralMixture instance with a single species.
    /// @param species The species that compose the mineral mixture
    explicit MineralMixture(const MineralSpecies& species);

    /// Destroy the MineralMixture instance.
    virtual ~MineralMixture();

    /// Calculate the state of the mineral mixture.
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    /// @param n The molar amounts of the species in the mixture (in units of mol)
    auto state(Temperature T, Pressure P, const Vector& n) const -> MineralMixtureState;
};

} // namespace Reaktoro
