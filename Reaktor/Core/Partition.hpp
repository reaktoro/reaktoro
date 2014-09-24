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

    /// Set the equilibrium species of the system using their indices
    /// @param indices The indices of the equilibrium species
    auto setEquilibriumSpecies(const Indices& indices) -> Partition&;

    /// Set the kinetic species of the system using their indices
    /// @param indices The indices of the kinetic species
    auto setKineticSpecies(const Indices& indices) -> Partition&;

    /// Set the inert species of the system using their indices
    /// @param indices The indices of the inert species
    auto setInertSpecies(const Indices& indices) -> Partition&;

    /// Get the indices of the equilibrium species
    auto indicesEquilibriumSpecies() const -> const Indices&;

    /// Get the indices of the kinetic species
    auto indicesKineticSpecies() const -> const Indices&;

    /// Get the indices of the inert species
    auto indicesInertSpecies() const -> const Indices&;

    /// Get the indices of the elements that compose the equilibrium species
    auto indicesEquilibriumElements() const -> const Indices&;

    /// Get the indices of the elements that compose the kinetic species
    auto indicesKineticElements() const -> const Indices&;

    /// Get the indices of the elements that compose the inert species
    auto indicesInertElements() const -> const Indices&;

//    /// Get the number of species in the system
//    auto numSpecies() const -> unsigned;
//
//    /// Get the number of equilibrium species in the system
//    auto numEquilibriumSpecies() const -> unsigned;
//
//    /// Get the number of kinetic species in the system
//    auto numKineticSpecies() const -> unsigned;
//
//    /// Get the number of inert species in the system
//    auto numInertSpecies() const -> unsigned;
//
//    /// Get the number of equilibrium elements in the system
//    auto numEquilibriumElements() const -> unsigned;
//
//    /// Get the number of kinetic elements in the system
//    auto numKineticElements() const -> unsigned;
//
//    /// Get the number of inert elements in the system
//    auto numInertElements() const -> unsigned;
//
//    /// Get the number of the phases that contains equilibrium species
//    auto numPhasesWithEquilibriumSpecies() const -> unsigned;
//
//    /// Get the number  of the phases that contains kinetic species
//    auto numPhasesWithKineticSpecies() const -> unsigned;
//
//    /// Get the number  of the phases that contains inert species
//    auto numPhasesWithInertSpecies() const -> unsigned;
//
//    /// Get the names of the species in the system
//    auto species() const -> const std::vector<std::string>&;
//
//    /// Get the names of the equilibrium species in the system
//    auto equilibriumSpecies() const -> const std::vector<std::string>&;
//
//    /// Get the names of the kinetic species in the system
//    auto kineticSpecies() const -> const std::vector<std::string>&;
//
//    /// Get the names of the inert species in the system
//    auto inertSpecies() const -> const std::vector<std::string>&;
//
//    /// Get the names of the equilibrium elements in the system
//    auto equilibriumElements() const -> const std::vector<std::string>&;
//
//    /// Get the names of the kinetic elements in the system
//    auto kineticElements() const -> const std::vector<std::string>&;
//
//    /// Get the names of the inert elements in the system
//    auto inertElements() const -> const std::vector<std::string>&;
//
//    /// Get the indices of the equilibrium species
//    auto idxEquilibriumSpecies() const -> const Indices&;
//
//    /// Get the indices of the kinetic species
//    auto idxKineticSpecies() const -> const Indices&;
//
//    /// Get the indices of the inert species
//    auto idxInertSpecies() const -> const Indices&;
//
//    /// Get the indices of the equilibrium elements
//    auto idxEquilibriumElements() const -> const Indices&;
//
//    /// Get the indices of the kinetic elements
//    auto idxKineticElements() const -> const Indices&;
//
//    /// Get the indices of the inert elements
//    auto idxInertElements() const -> const Indices&;
//
//    /// Get the index of a species in the system
//    auto idxSpecies(const std::string& species) const -> Index;
//
//    /// Get the indices of a list of species in the system
//    auto idxSpecies(const std::vector<std::string>& species) const -> Indices;
//
//    /// Get the index of an equilibrium species in the system
//    auto idxEquilibriumSpecies(const std::string& species) const -> Index;
//
//    /// Get the indices of a list of equilibrium species in the system
//    auto idxEquilibriumSpecies(const std::vector<std::string>& species) const -> Indices;
//
//    /// Get the index of a kinetic species in the system
//    auto idxKineticSpecies(const std::string& species) const -> Index;
//
//    /// Get the indices of a list of kinetic species in the system
//    auto idxKineticSpecies(const std::vector<std::string>& species) const -> Indices;
//
//    /// Get the index of a inert species in the system
//    auto idxInertSpecies(const std::string& species) const -> Index;
//
//    /// Get the indices of a list of inert species in the system
//    auto idxInertSpecies(const std::vector<std::string>& species) const -> Indices;
//
//    /// Get the index of an equilibrium element in the system
//    auto idxEquilibriumElement(const std::string& species) const -> Index;
//
//    /// Get the indices of a list of equilibrium elements in the system
//    auto idxEquilibriumElements(const std::vector<std::string>& species) const -> Indices;
//
//    /// Get the index of a kinetic element in the system
//    auto idxKineticElement(const std::string& species) const -> Index;
//
//    /// Get the indices of a list of kinetic elements in the system
//    auto idxKineticElements(const std::vector<std::string>& species) const -> Indices;
//
//    /// Get the index of a inert element in the system
//    auto idxInertElement(const std::string& species) const -> Index;
//
//    /// Get the indices of a list of inert elements in the system
//    auto idxInertElements(const std::vector<std::string>& species) const -> Indices;
//
//    /// Get the indices of the phases that contains equilibrium species
//    auto idxPhasesWithEquilibriumSpecies() const -> const Indices&;
//
//    /// Get the indices of the phases that contains kinetic species
//    auto idxPhasesWithKineticSpecies() const -> const Indices&;
//
//    /// Get the indices of the phases that contains inert species
//    auto idxPhasesWithInertSpecies() const -> const Indices&;
//
//    /// Set the rows of a vector whose indices correspond to the equilibrium species
//    ///
//    /// This method is useful to set the composition of the equilibrium species.
//    ///
//    /// **Usage**
//    ///
//    /// In the example below, let `n` denote the compositional vector of the system,
//    /// and `ne` the compositional vector of the equilibrium species.
//    /// Updating the composition of the equilibrium species in `n` can then be done
//    /// with the call:
//    /// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    /// Partition.setEquilibriumRows(ne, n);
//    /// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    ///
//    /// @param values The values to be used when setting the equilibrium rows of `vec`
//    /// @param vec The vector whose equilibrium rows will be set
//    ///
//    /// @see setKineticRows, setInertRows
//    auto setEquilibriumRows(const Vector& values, Vector& vec) const -> void;
//
//    /// Set the rows of a vector whose indices correspond to the kinetic species
//    ///
//    /// This method is useful to set the composition of the kinetic species.
//    ///
//    /// **Usage**
//    ///
//    /// In the example below, let `n` denote the compositional vector of the system,
//    /// and `nk` the compositional vector of the kinetic species.
//    /// Updating the composition of the kinetic species in `n` can then be done
//    /// with the call:
//    /// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    /// Partition.setEquilibriumRows(nk, n);
//    /// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    ///
//    /// @param values The values to be used when setting the kinetic rows of `vec`
//    /// @param vec The vector whose kinetic rows will be set
//    ///
//    /// @see setEquilibriumRows, setInertRows
//    auto setKineticRows(const Vector& values, Vector& vec) const -> void;
//
//    /// Set the rows of a vector whose indices correspond to the inert species
//    ///
//    /// This method is useful to set the composition of the inert species.
//    ///
//    /// **Usage**
//    ///
//    /// In the example below, let `n` denote the compositional vector of the system,
//    /// and `ni` the compositional vector of the inert species.
//    /// Updating the composition of the inert species in `n` can then be done
//    /// with the call:
//    /// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    /// Partition.setEquilibriumRows(ni, n);
//    /// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    ///
//    /// @param values The values to be used when setting the inert rows of `vec`
//    /// @param vec The vector whose inert rows will be set
//    ///
//    /// @see setEquilibriumRows, setKineticRows
//    auto setInertRows(const Vector& values, Vector& vec) const -> void;
//
//    /// Get the rows of a vector whose indices correspond to the equilibrium species
//    /// @param vec The vector whose equilibrium rows will be extracted
//    auto equilibriumRows(const Vector& vec) const -> Vector;
//
//    /// Get the rows of a vector whose indices correspond to the kinetic species
//    /// @param vec The vector whose kinetic rows will be extracted
//    auto kineticRows(const Vector& vec) const -> Vector;
//
//    /// Get the rows of a vector whose indices correspond to the inert species
//    /// @param vec The vector whose inert rows will be extracted
//    auto inertRows(const Vector& vec) const -> Vector;
//
//    /// Get the columns of a matrix whose indices correspond to the equilibrium species
//    /// @param vec The matrix whose equilibrium columns will be extracted
//    auto equilibriumCols(const Matrix& mat) const -> Matrix;
//
//    /// Get the columns of a matrix whose indices correspond to the kinetic species
//    /// @param vec The matrix whose kinetic columns will be extracted
//    auto kineticCols(const Matrix& mat) const -> Matrix;
//
//    /// Get the columns of a matrix whose indices correspond to the inert species
//    /// @param vec The matrix whose inert columns will be extracted
//    auto inertCols(const Matrix& mat) const -> Matrix;
//
//    /// Get the rows and columns of a matrix whose indices correspond to the equilibrium species
//    /// @param vec The matrix whose equilibrium rows and columns will be extracted
//    auto equilibriumRowsCols(const Matrix& mat) const -> Matrix;
//
//    /// Get the rows and columns of a matrix whose indices correspond to the kinetic species
//    /// @param vec The matrix whose kinetic rows and columns will be extracted
//    auto kineticRowsCols(const Matrix& mat) const -> Matrix;
//
//    /// Get the rows and columns of a matrix whose indices correspond to the inert species
//    /// @param vec The matrix whose inert rows and columns will be extracted
//    auto inertRowsCols(const Matrix& mat) const -> Matrix;
//
//    /// Get the columns of the formula matrix whose indices correspond to the equilibrium species
//    /// @param system The chemical system instance
//    /// @see Chemical, ChemicalSystem::formulaMatrix
//    auto equilibriumFormulaMatrix(const ChemicalSystem& system) const -> Matrix;
//
//    /// Get the columns of the formula matrix whose indices correspond to the kinetic species
//    /// @param system The chemical system instance
//    /// @see ChemicalSystem, ChemicalSystem::formulaMatrix
//    auto kineticFormulaMatrix(const ChemicalSystem& system) const -> Matrix;
//
//    /// Get the columns of the formula matrix whose indices correspond to the inert species
//    /// @param system The chemical system instance
//    /// @see ChemicalSystem, ChemicalSystem::formulaMatrix
//    auto inertFormulaMatrix(const ChemicalSystem& system) const -> Matrix;
//
//    /// Get the columns of the stoichiometric matrix whose indices correspond to the equilibrium species
//    /// @param reactions The reaction system instance
//    /// @see ReactionSystem, ReactionSystem::stoichiometricMatrix
//    auto equilibriumStoichiometricMatrix(const ReactionSystem& reactions) const -> Matrix;
//
//    /// Get the columns of the stoichiometric matrix whose indices correspond to the kinetic species
//    /// @param reactions The reaction system instance
//    /// @see ReactionSystem, ReactionSystem::stoichiometricMatrix
//    auto kineticStoichiometricMatrix(const ReactionSystem& reactions) const -> Matrix;
//
//    /// Get the columns of the stoichiometric matrix whose indices correspond to the inert species
//    /// @param reactions The reaction system instance
//    /// @see ReactionSystem, ReactionSystem::stoichiometricMatrix
//    auto inertStoichiometricMatrix(const ReactionSystem& reactions) const -> Matrix;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktor
