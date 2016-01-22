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

// C++ includes
#include <memory>

// Reaktoro includes
#include <Reaktoro/Core/Phase.hpp>

namespace Reaktoro {

// Forward declarations
class MineralMixture;
class MineralSpecies;

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

    /// Return the MineralMixture instance
    auto mixture() const -> const MineralMixture&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
