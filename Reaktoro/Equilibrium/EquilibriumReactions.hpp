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
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Math/LU.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;
class Partition;
class ReactionEquation;

/// A class that generates a system of equilibrium reactions written in terms of master and secondary species.
class EquilibriumReactions
{
public:
    /// Construct an EquilibriumReactions instance.
    EquilibriumReactions(const ChemicalSystem& system);

    /// Construct an EquilibriumReactions instance
    EquilibriumReactions(const ChemicalSystem& system, const Partition& partition);

    /// Construct a copy of an EquilibriumReactions instance
    EquilibriumReactions(const EquilibriumReactions& other);

    /// Destroy this EquilibriumReactions instance.
    virtual ~EquilibriumReactions();

    /// Assign other EquilibriumReactions instance to this.
    auto operator=(EquilibriumReactions other) -> EquilibriumReactions&;

    /// Return the chemical system for which the equilibrium reactions were defined.
    auto system() const -> const ChemicalSystem&;

    /// Return the partition of the chemical system for which the equilibrium reactions were defined.
    auto partition() const -> const Partition&;

    /// Set the master species manually.
    /// @param ispecies The global indices of the master species.
    auto setMasterSpecies(Indices ispecies) -> void;

    /// Set the master species manually.
    /// @param species The names of the master species.
    auto setMasterSpecies(std::vector<std::string> species) -> void;

    /// Return the indices of the master species.
    /// The master species are those that serve as building blocks for the secondary species.
    auto indicesMasterSpecies() const -> Indices;

    /// Return the indices of the secondary species.
    /// The secondary species are those that are constructed from master species.
    auto indicesSecondarySpecies() const -> Indices;

    /// Return the equations of the equilibrium reactions.
    auto equations() const -> std::vector<ReactionEquation>;

    /// Return the stoichiometric matrix of the reactions.
    auto stoichiometricMatrix() const -> Matrix;

    /// Return the LU decomposition of the formula matrix `A`.
    auto lu() const -> const LU&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
