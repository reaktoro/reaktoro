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
#include <string>

// Reaktoro includes
#include <Reaktoro/Equilibrium/EquilibriumState.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;

/// Provides a computational representation of the kinetic state of a multiphase chemical system.
/// @see EquilibriumState, ChemicalState, ChemicalSystem
/// @ingroup Core
class KineticState : public EquilibriumState
{
public:
    /// Construct a default KineticState instance
    KineticState();

    /// Construct an KineticState instance using a ChemicalSystem instance.
    /// @param system The chemical system instance
    explicit KineticState(const ChemicalSystem& system);

    /// Construct an KineticState instance using a ChemicalState instance.
    /// @param state The chemical state instance
    explicit KineticState(const ChemicalState& state);

    /// Construct an KineticState instance using a EquilibriumState instance.
    /// @param state The chemical state instance
    KineticState(const EquilibriumState& state);

    /// Construct a copy of an KineticState instance
    KineticState(const KineticState& other);

    /// Destroy the instance
    virtual ~KineticState();

    /// Assign an KineticState instance to this instance
    auto operator=(KineticState other) -> KineticState&;

    using EquilibriumState::operator=;

    /// Output the KineticState instance to a file.
    auto output(std::string filename) -> void;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

/// Outputs an KineticState instance.
auto operator<<(std::ostream& out, const KineticState& state) -> std::ostream&;

} // namespace Reaktoro
