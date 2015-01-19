// Reaktor is a C++ library for computational reaction modelling.
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
#include <memory>
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/Matrix.hpp>

namespace Reaktor {

// Forward declarations
class ChemicalSystem;

/// Provide a computational representation of the Partition of a chemical system
///
/// A chemical system can be partitioned into *equilibrium*, *kinetic* and *inert species*.
///
/// The equilibrium species are the species whose composition is governed by chemical
/// equilibrium. In other words, their composition is calculated by the minimization of
/// their Gibbs energy subject to some equilibrium constraints (e.g., mass-balance
/// constraints).
///
/// The kinetic species are the species whose composition is governed by chemical
/// kinetics. By solving a system of ordinary differential equations that model the
/// kinetics of a system of reactions, the composition of the kinetic species can be
/// traced with time. The composition of the equilibrium species with time is calculated
/// with chemical equilibrium calculations with equilibrium constraints accounting for
/// the kinetic variation of the molar abundance of the chemical elements in the
/// equilibrium Partition.
///
/// The inert species are the species whose composition is invariable.
///
/// @see ChemicalSystem
/// @ingroup Core
class Partition
{
public:
    /// Construct a Partition instance
    /// @param system The chemical system instance
    /// @see ChemicalSystem
    Partition(const ChemicalSystem& system);

    /// Construct a copy of a Partition instance
    Partition(const Partition& other);

    /// Destroy the instance
    virtual ~Partition();

    /// Assign a Partition instance to this instance
    auto operator=(Partition other) -> Partition&;

    /// Set the equilibrium species of the chemical system
    auto setEquilibriumSpecies(const Indices& ispecies) -> void;

    /// Set the equilibrium species of the chemical system
    auto setEquilibriumSpecies(const std::vector<std::string>& species) -> void;

    /// Set the kinetic species of the chemical system
    auto setKineticSpecies(const Indices& species) -> void;

    /// Set the kinetic species of the chemical system
    auto setKineticSpecies(const std::vector<std::string>& species) -> void;

    /// Set the inert species of the chemical system
    auto setInertSpecies(const Indices& species) -> void;

    /// Set the inert species of the chemical system
    auto setInertSpecies(const std::vector<std::string>& species) -> void;

    /// Get the indices of the equilibrium species in the partition
    auto indicesEquilibriumSpecies() const -> const Indices&;

    /// Get the indices of the kinetic species in the partition
    auto indicesKineticSpecies() const -> const Indices&;

    /// Get the indices of the inert species in the partition
    auto indicesInertSpecies() const -> const Indices&;

    /// Get the indices of the elements in the equilibrium partition
    auto indicesElementsInEquilibriumSpecies() const -> const Indices&;

    /// Get the indices of the elements in the kinetic partition
    auto indicesElementsInKineticSpecies() const -> const Indices&;

    /// Get the indices of the elements in the inert partition
    auto indicesElementsInInertSpecies() const -> const Indices&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktor
