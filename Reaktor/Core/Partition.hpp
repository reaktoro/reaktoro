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
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

// Forward declarations
class Multiphase;

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
/// @see Multiphase
/// @ingroup Core
class Partition
{
public:
    /// Construct a default Partition instance
    Partition() = delete;

    /// Construct a Partition instance
    /// @param multiphase The multiphase system
    /// @see Multiphase
    explicit Partition(const Multiphase& multiphase);

    /// Construct a copy of a Partition instance
    Partition(const Partition& other);

    /// Destroy the instance
    virtual ~Partition();

    /// Assign a Partition instance to this instance
    auto operator=(Partition other) -> Partition&;

    /// Set the indices of the equilibrium species in the partition
    /// @param indices The indices of the equilibrium species
    auto setIndicesEquilibriumSpecies(const Indices& indices) -> Partition&;

    /// Set the indices of the kinetic species in the partition
    /// @param indices The indices of the kinetic species
    auto setIndicesKineticSpecies(const Indices& indices) -> Partition&;

    /// Set the indices of the inert species in the partition
    /// @param indices The indices of the inert species
    auto setIndicesInertSpecies(const Indices& indices) -> Partition&;

    /// Get the indices of the equilibrium species in the partition
    auto indicesEquilibriumSpecies() const -> const Indices&;

    /// Get the indices of the kinetic species in the partition
    auto indicesKineticSpecies() const -> const Indices&;

    /// Get the indices of the inert species in the partition
    auto indicesInertSpecies() const -> const Indices&;

    /// Get the indices of the elements that compose the equilibrium species in the partition
    auto indicesEquilibriumElements() const -> const Indices&;

    /// Get the indices of the elements that compose the kinetic species in the partition
    auto indicesKineticElements() const -> const Indices&;

    /// Get the indices of the elements that compose the inert species in the partition
    auto indicesInertElements() const -> const Indices&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktor
