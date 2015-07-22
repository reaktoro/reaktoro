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
#include <memory>
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;

/// Provide a computational representation of the Partition of a chemical system.
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
    /// Construct a default Partition instance
    Partition();

    /// Construct a Partition instance
    /// @param system The chemical system instance
    /// @see ChemicalSystem
    Partition(const ChemicalSystem& system);

    /// Construct a Partition instance using a formatted string
    Partition(const ChemicalSystem& system, std::string partition);

    /// Construct a copy of a Partition instance
    Partition(const Partition& other);

    /// Destroy the instance
    virtual ~Partition();

    /// Assign a Partition instance to this instance
    auto operator=(Partition other) -> Partition&;

    /// Set the partition of the chemical system using a formatted string
    auto set(std::string partition) -> void;

    /// Set the equilibrium species of the chemical system
    auto setEquilibriumSpecies(const Indices& ispecies) -> void;

    /// Set the equilibrium species of the chemical system
    auto setEquilibriumSpecies(const std::vector<std::string>& species) -> void;

    /// Set the equilibrium species of the chemical system as the species in given phases
    auto setEquilibriumPhases(const Indices& iphases) -> void;

    /// Set the equilibrium species of the chemical system as the species in given phases
    auto setEquilibriumPhases(const std::vector<std::string>& phases) -> void;

    /// Set the kinetic species of the chemical system
    auto setKineticSpecies(const Indices& species) -> void;

    /// Set the kinetic species of the chemical system
    auto setKineticSpecies(const std::vector<std::string>& species) -> void;

    /// Set the kinetic species of the chemical system as the species in given phases
    auto setKineticPhases(const Indices& iphases) -> void;

    /// Set the kinetic species of the chemical system as the species in given phases
    auto setKineticPhases(const std::vector<std::string>& phases) -> void;

    /// Set the inert species of the chemical system
    auto setInertSpecies(const Indices& species) -> void;

    /// Set the inert species of the chemical system
    auto setInertSpecies(const std::vector<std::string>& species) -> void;

    /// Set the inert species of the chemical system as the species in given phases
    auto setInertPhases(const Indices& iphases) -> void;

    /// Set the inert species of the chemical system as the species in given phases
    auto setInertPhases(const std::vector<std::string>& phases) -> void;

    /// Return the number of equilibrium species in the partition
    auto numEquilibriumSpecies() const -> unsigned;

    /// Return the number of kinetic species in the partition
    auto numKineticSpecies() const -> unsigned;

    /// Return the number of inert species in the partition
    auto numInertSpecies() const -> unsigned;

    /// Return the number of equilibrium elements in the partition
    auto numEquilibriumElements() const -> unsigned;

    /// Return the number of kinetic elements in the partition
    auto numKineticElements() const -> unsigned;

    /// Return the number of inert elements in the partition
    auto numInertElements() const -> unsigned;

    /// Return the indices of the equilibrium species in the partition
    auto indicesEquilibriumSpecies() const -> const Indices&;

    /// Return the indices of the kinetic species in the partition
    auto indicesKineticSpecies() const -> const Indices&;

    /// Return the indices of the inert species in the partition
    auto indicesInertSpecies() const -> const Indices&;

    /// Return the indices of the elements in the equilibrium partition
    auto indicesEquilibriumElements() const -> const Indices&;

    /// Return the indices of the elements in the kinetic partition
    auto indicesKineticElements() const -> const Indices&;

    /// Return the indices of the elements in the inert partition
    auto indicesInertElements() const -> const Indices&;

    /// Return the formula matrix of the equilibrium species with respect to the equilibrium elements
    auto formulaMatrixEquilibriumSpecies() const -> const Matrix&;

    /// Return the formula matrix of the kinetic species with respect to the kinetic elements
    auto formulaMatrixKineticSpecies() const -> const Matrix&;

    /// Return the formula matrix of the inert species with respect to the inert elements
    auto formulaMatrixInertSpecies() const -> const Matrix&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
