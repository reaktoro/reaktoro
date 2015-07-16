// Reaktoro is a C++ library for computational reaction modelling.
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

#pragma once

// C++ includes
#include <string>
#include <memory>

// Reaktoro includes
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Thermodynamics/Activity/MineralActivity.hpp>

namespace Reaktoro {

// Forward declarations
class MineralMixture;

/// Class that defines an mineral phase
class MineralPhase : public Phase
{
public:
    /// Construct a default MineralPhase instance.
    MineralPhase();

    /// Construct a copy of a MineralPhase instance
    MineralPhase(const MineralPhase& other);

    /// Construct a MineralPhase instance with given mineral mixture.
    explicit MineralPhase(const MineralMixture& mixture);

    /// Construct a MineralPhase instance with given species.
    explicit MineralPhase(const MineralSpecies& species);

    /// Destroy the MineralPhase instance.
    virtual ~MineralPhase();

    /// Assign a MineralPhase instance to this
    auto operator=(MineralPhase other) -> MineralPhase&;

    /// Set the chemical model of the phase with the ideal solution model.
    auto setChemicalModelIdeal() -> void;

    /// Set the activity model of a species.
    /// @param species The name of the species
    /// @param activity The activity function
    /// @see MineralActivity
    auto setActivityModel(const std::string& species, const MineralActivityFunction& activity) -> void;

    /// Set the activity model of the species to be the ideal one.
    /// @param species The name of species to have its activity model set
    auto setActivityModelIdeal(const std::string& species) -> void;

    /// Return the MineralMixture instance
    auto mixture() const -> const MineralMixture&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
