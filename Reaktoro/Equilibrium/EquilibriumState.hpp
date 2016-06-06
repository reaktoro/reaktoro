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
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;

/// Provides a computational representation of the equilibrium state of a multiphase chemical system.
/// @see ChemicalState, ChemicalSystem
/// @ingroup Core
class EquilibriumState : public ChemicalState
{
public:
    /// Construct a default EquilibriumState instance
    EquilibriumState();

    /// Construct an EquilibriumState instance using a ChemicalSystem instance.
    /// @param system The chemical system instance
    explicit EquilibriumState(const ChemicalSystem& system);

    /// Construct an EquilibriumState instance using a ChemicalState instance.
    /// @param state The chemical state instance
    EquilibriumState(const ChemicalState& state);

    /// Construct a copy of an EquilibriumState instance
    EquilibriumState(const EquilibriumState& other);

    /// Destroy the instance
    virtual ~EquilibriumState();

    /// Assign an EquilibriumState instance to this instance
    auto operator=(EquilibriumState other) -> EquilibriumState&;

    using ChemicalState::operator=;

    /// Set the dual potentials of the species (in units of J/mol)
    /// The dual potentials of the species are the Lagrange multipliers with
    /// respect to the positive bound constraints on the molar amounts of the
    /// species in a chemical equilibrium calculation. They can be seen as
    /// measures of stability of a species at equilibrium, with values closer
    /// to zero meaning more stability.
    /// @param values The Lagrange multipliers with respect to the positive constraints.
    auto setSpeciesDualPotentials(const Vector& values) -> void;

    /// Set the dual potentials of the elements (in units of J/mol)
    /// The dual potentials of the elements are the Lagrange multipliers with
    /// respect to the balance constraints on the molar amounts of the elements.
    /// They can be seen as dual chemical potential of elements.
    /// @param values The Lagrange multipliers with respect to the balance constraints.
    auto setElementDualPotentials(const Vector& values) -> void;

    /// Return the dual potentials of the species (in units of J/mol)
    auto speciesDualPotentials() const -> const Vector&;

    /// Return the dual potentials of the elements (in units of J/mol)
    auto elementDualPotentials() const -> const Vector&;

    /// Return the stability indices of the phases with respect to chemical equilibrium.
    /// The stability index of a stable phase at chemical equilibrium should
    /// be zero or very close to zero. A negative stability index indicates
    /// that the corresponding phase is under-saturated, while a positive index
    /// indicates the phase is over-saturated.
    auto phaseStabilityIndices() const -> Vector;

    /// Output the EquilibriumState instance to a file.
    auto output(std::string filename) -> void;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

/// Outputs an EquilibriumState instance.
auto operator<<(std::ostream& out, const EquilibriumState& state) -> std::ostream&;

} // namespace Reaktoro
