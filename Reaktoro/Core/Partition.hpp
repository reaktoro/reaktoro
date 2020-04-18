// // Reaktoro is a unified framework for modeling chemically reactive systems.
// //
// // Copyright (C) 2014-2020 Allan Leal
// //
// // This library is free software; you can redistribute it and/or
// // modify it under the terms of the GNU Lesser General Public
// // License as published by the Free Software Foundation; either
// // version 2.1 of the License, or (at your option) any later version.
// //
// // This library is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// // Lesser General Public License for more details.
// //
// // You should have received a copy of the GNU Lesser General Public License
// // along with this library. If not, see <http://www.gnu.org/licenses/>.

// #pragma once

// // C++ includes
// #include <memory>
// #include <string>
// #include <vector>

// // Reaktoro includes
// #include <Reaktoro/Common/Index.hpp>
// #include <Reaktoro/Math/Matrix.hpp>

// namespace Reaktoro {

// // Forward declarations
// class ChemicalSystem;

// /// Provide a computational representation of the Partition of a chemical system.
// ///
// /// A chemical system can be partitioned into *equilibrium*, *kinetic* and *inert species*.
// ///
// /// The equilibrium species are the species whose composition is governed by chemical
// /// equilibrium. In other words, their composition is calculated by the minimization of
// /// their Gibbs energy subject to some equilibrium constraints (e.g., mass-balance
// /// constraints).
// ///
// /// The kinetic species are the species whose composition is governed by chemical
// /// kinetics. By solving a system of ordinary differential equations that model the
// /// kinetics of a system of reactions, the composition of the kinetic species can be
// /// traced with time. The composition of the equilibrium species with time is calculated
// /// with chemical equilibrium calculations with equilibrium constraints accounting for
// /// the kinetic variation of the molar abundance of the chemical elements in the
// /// equilibrium Partition.
// ///
// /// The inert species are the species whose composition is invariable.
// ///
// /// @see ChemicalSystem
// /// @ingroup Core
// class Partition
// {
// public:
//     /// Construct a default Partition instance
//     Partition();

//     /// Construct a Partition instance
//     /// @param system The chemical system instance
//     /// @see ChemicalSystem
//     Partition(const ChemicalSystem& system);

//     /// Set the equilibrium species of the chemical system
//     auto setEquilibriumSpecies(const Indices& ispecies) -> void;

//     /// Set the equilibrium species of the chemical system
//     auto setEquilibriumSpecies(const std::vector<std::string>& species) -> void;

//     /// Set the equilibrium species of the chemical system as the species in given phases
//     auto setEquilibriumPhases(const Indices& iphases) -> void;

//     /// Set the equilibrium species of the chemical system as the species in given phases
//     auto setEquilibriumPhases(const std::vector<std::string>& phases) -> void;

//     /// Set the kinetic species of the chemical system
//     auto setKineticSpecies(const Indices& species) -> void;

//     /// Set the kinetic species of the chemical system
//     auto setKineticSpecies(const std::vector<std::string>& species) -> void;

//     /// Set the kinetic species of the chemical system as the species in given phases
//     auto setKineticPhases(const Indices& iphases) -> void;

//     /// Set the kinetic species of the chemical system as the species in given phases
//     auto setKineticPhases(const std::vector<std::string>& phases) -> void;

//     /// Set the inert species of the chemical system
//     auto setInertSpecies(const Indices& species) -> void;

//     /// Set the inert species of the chemical system
//     auto setInertSpecies(const std::vector<std::string>& species) -> void;

//     /// Set the inert species of the chemical system as the species in given phases
//     auto setInertPhases(const Indices& iphases) -> void;

//     /// Set the inert species of the chemical system as the species in given phases
//     auto setInertPhases(const std::vector<std::string>& phases) -> void;

//     /// Set the fluid phases of the chemical system.
//     auto setFluidPhases(const Indices& indices) -> void;

//     /// Set the fluid phases of the chemical system.
//     auto setFluidPhases(const std::vector<std::string>& names) -> void;

//     /// Set the solid phases of the chemical system.
//     auto setSolidPhases(const Indices& indices) -> void;

//     /// Set the solid phases of the chemical system.
//     auto setSolidPhases(const std::vector<std::string>& names) -> void;

//     /// Return the chemical system.
//     auto system() const -> const ChemicalSystem&;

//     /// Return the number of phases in the fluid partition.
//     auto numFluidPhases() const -> unsigned;

//     /// Return the number of species in the fluid partition.
//     auto numFluidSpecies() const -> unsigned;

//     /// Return the number of phases in the solid partition.
//     auto numSolidPhases() const -> unsigned;

//     /// Return the number of species in the solid partition.
//     auto numSolidSpecies() const -> unsigned;

//     /// Return the number of species in the equilibrium partition.
//     auto numEquilibriumSpecies() const -> unsigned;

//     /// Return the number of species in the equilibrium-fluid partition.
//     auto numEquilibriumFluidSpecies() const -> unsigned;

//     /// Return the number of species in the equilibrium-solid partition.
//     auto numEquilibriumSolidSpecies() const -> unsigned;

//     /// Return the number of species in the kinetic partition.
//     auto numKineticSpecies() const -> unsigned;

//     /// Return the number of species in the kinetic-fluid partition.
//     auto numKineticFluidSpecies() const -> unsigned;

//     /// Return the number of species in the kinetic-solid partition.
//     auto numKineticSolidSpecies() const -> unsigned;

//     /// Return the number of species in the inert partition.
//     auto numInertSpecies() const -> unsigned;

//     /// Return the number of species in the inert-fluid partition.
//     auto numInertFluidSpecies() const -> unsigned;

//     /// Return the number of species in the inert-solid partition.
//     auto numInertSolidSpecies() const -> unsigned;

//     /// Return the number of elements in the equilibrium partition.
//     auto numEquilibriumElements() const -> unsigned;

//     /// Return the number of elements in the equilibrium-fluid partition.
//     auto numEquilibriumFluidElements() const -> unsigned;

//     /// Return the number of elements in the equilibrium-solid partition.
//     auto numEquilibriumSolidElements() const -> unsigned;

//     /// Return the number of elements in the kinetic partition.
//     auto numKineticElements() const -> unsigned;

//     /// Return the number of elements in the kinetic-fluid partition.
//     auto numKineticFluidElements() const -> unsigned;

//     /// Return the number of elements in the kinetic-solid partition.
//     auto numKineticSolidElements() const -> unsigned;

//     /// Return the number of elements in the inert partition.
//     auto numInertElements() const -> unsigned;

//     /// Return the number of elements in the inert-fluid partition.
//     auto numInertFluidElements() const -> unsigned;

//     /// Return the number of elements in the inert-solid partition.
//     auto numInertSolidElements() const -> unsigned;

//     /// Return the indices of the phases in the fluid partition.
//     auto indicesFluidPhases() const -> const Indices&;

//     /// Return the indices of the species in the fluid partition.
//     auto indicesFluidSpecies() const -> const Indices&;

//     /// Return the indices of the phases in the solid partition.
//     auto indicesSolidPhases() const -> const Indices&;

//     /// Return the indices of the species in the solid partition.
//     auto indicesSolidSpecies() const -> const Indices&;

//     /// Return the indices of the species in the equilibrium partition.
//     auto indicesEquilibriumSpecies() const -> const Indices&;

//     /// Return the indices of the species in the equilibrium-fluid partition.
//     auto indicesEquilibriumFluidSpecies() const -> const Indices&;

//     /// Return the indices of the species in the equilibrium-solid partition.
//     auto indicesEquilibriumSolidSpecies() const -> const Indices&;

//     /// Return the indices of the species in the kinetic partition.
//     auto indicesKineticSpecies() const -> const Indices&;

//     /// Return the indices of the species in the kinetic-fluid partition.
//     auto indicesKineticFluidSpecies() const -> const Indices&;

//     /// Return the indices of the species in the kinetic-solid partition.
//     auto indicesKineticSolidSpecies() const -> const Indices&;

//     /// Return the indices of the species in the inert partition.
//     auto indicesInertSpecies() const -> const Indices&;

//     /// Return the indices of the species in the inert-fluid partition.
//     auto indicesInertFluidSpecies() const -> const Indices&;

//     /// Return the indices of the species in the inert-solid partition.
//     auto indicesInertSolidSpecies() const -> const Indices&;

//     /// Return the indices of the elements in the equilibrium partition.
//     auto indicesEquilibriumElements() const -> const Indices&;

//     /// Return the indices of the elements in the equilibrium-fluid partition.
//     auto indicesEquilibriumFluidElements() const -> const Indices&;

//     /// Return the indices of the elements in the equilibrium-solid partition.
//     auto indicesEquilibriumSolidElements() const -> const Indices&;

//     /// Return the indices of the elements in the kinetic partition.
//     auto indicesKineticElements() const -> const Indices&;

//     /// Return the indices of the elements in the kinetic-fluid partition.
//     auto indicesKineticFluidElements() const -> const Indices&;

//     /// Return the indices of the elements in the kinetic-solid partition.
//     auto indicesKineticSolidElements() const -> const Indices&;

//     /// Return the indices of the elements in the inert partition.
//     auto indicesInertElements() const -> const Indices&;

//     /// Return the indices of the elements in the inert-fluid partition.
//     auto indicesInertFluidElements() const -> const Indices&;

//     /// Return the indices of the elements in the inert-solid partition.
//     auto indicesInertSolidElements() const -> const Indices&;

//     /// Return the formula matrix of the equilibrium partition.
//     auto formulaMatrixEquilibriumPartition() const -> MatrixXdConstRef;

//     /// Return the formula matrix of the equilibrium-fluid partition.
//     auto formulaMatrixEquilibriumFluidPartition() const -> MatrixXdConstRef;

//     /// Return the formula matrix of the equilibrium-solid partition.
//     auto formulaMatrixEquilibriumSolidPartition() const -> MatrixXdConstRef;

//     /// Return the formula matrix of the kinetic partition.
//     auto formulaMatrixKineticPartition() const -> MatrixXdConstRef;

//     /// Return the formula matrix of the kinetic-fluid partition.
//     auto formulaMatrixKineticFluidPartition() const -> MatrixXdConstRef;

//     /// Return the formula matrix of the kinetic-solid partition.
//     auto formulaMatrixKineticSolidPartition() const -> MatrixXdConstRef;

//     /// Return the formula matrix of the inert partition.
//     auto formulaMatrixInertPartition() const -> MatrixXdConstRef;

//     /// Return the formula matrix of the inert-fluid partition.
//     auto formulaMatrixInertFluidPartition() const -> MatrixXdConstRef;

//     /// Return the formula matrix of the inert-solid partition.
//     auto formulaMatrixInertSolidPartition() const -> MatrixXdConstRef;

// private:
//     struct Impl;

//     std::shared_ptr<Impl> pimpl;
// };

// } // namespace Reaktoro
