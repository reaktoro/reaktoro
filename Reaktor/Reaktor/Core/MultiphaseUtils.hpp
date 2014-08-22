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
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/Vector.hpp>
#include <Reaktor/Common/Matrix.hpp>

namespace Reaktor {

// Forward declarations
class Multiphase;
class ScalarResult;
class VectorResult;

/// Get the number of elements in a multiphase system
auto numElements(const Multiphase& multiphase) -> unsigned;

/// Get the number of species in a multiphase system
auto numSpecies(const Multiphase& multiphase) -> unsigned;

/// Get the number of phases in a multiphase system
auto numPhases(const Multiphase& multiphase) -> unsigned;

/// Determine if a multiphase system contains an element
/// @param system The multiphase system
/// @param species The name of the element
auto containsElement(const Multiphase& multiphase, const std::string& element) -> bool;

/// Determine if a multiphase system contains a species
/// @param system The multiphase system
/// @param species The name of the species
auto containsSpecies(const Multiphase& multiphase, const std::string& species) -> bool;

/// Determine if a multiphase system contains a phase
/// @param system The multiphase system
/// @param species The name of the phase
auto containsPhase(const Multiphase& multiphase, const std::string& phase) -> bool;

/// Get the index of an element in a multiphase system
/// @param system The multiphase system
/// @param name The name of the element
/// @return The index of the element if found. Otherwise the number of elements in the chemical system
auto indexElement(const Multiphase& multiphase, const std::string& name) -> Index;

/// Get the indices of some elements in a multiphase system
/// @param system The multiphase system
/// @param names The names of the elements
/// @return The indices of the elements if found, or the number of elements in the system otherwise
auto indicesElements(const Multiphase& multiphase, const std::vector<std::string>& names) -> Indices;

/// Get the index of a species in a multiphase system
/// @param system The multiphase system
/// @param name The name of the species
/// @return The index of the species if found, or the number of species in the system otherwise
auto indexSpecies(const Multiphase& multiphase, const std::string& name) -> Index;

/// Get the indices of some species in a multiphase system
/// @param system The multiphase system
/// @param names The names of the species
/// @return The indices of the species if found, or the number of species in the system otherwise
auto indicesSpecies(const Multiphase& multiphase, const std::vector<std::string>& names) -> Indices;

/// Get the index of a phase in a multiphase system
/// @param system The multiphase system
/// @param name The name of the phase
/// @return The index of the phase if found. Otherwise the number of phases in the chemical system
auto indexPhase(const Multiphase& multiphase, const std::string& name) -> Index;

/// Get the indices of some phases in a multiphase system
/// @param system The multiphase system
/// @param names The names of the phases
/// @return The indices of the phases (the number of phases in the chemical system for the not found phases)
auto indicesPhases(const Multiphase& multiphase, const std::vector<std::string>& names) -> Indices;

/// Get the index of the first species of a given phase in a multiphase system
/// @param system The multiphase system
/// @param iphase The index of the phase in the chemical system
auto indexFirstSpeciesInPhase(const Multiphase& multiphase, const Index& iphase) -> Index;

/// Get the index of the last species of a given phase in a multiphase system
/// @param system The multiphase system
/// @param iphase The index of the phase in the chemical system
auto indexLastSpeciesInPhase(const Multiphase& multiphase, const Index& iphase) -> Index;

/// Get the indices of the elements that compose a given species in a multiphase system
/// @param system The multiphase system
/// @param ispecies The index of the species in the chemical system
auto indicesElementsInSpecies(const Multiphase& multiphase, const Index& ispecies) -> Indices;

/// Get the indices of the species that compose a given phase in a multiphase system
/// @param system The multiphase system
/// @param iphase The index of the phase in the chemical system
auto indicesSpeciesInPhase(const Multiphase& multiphase, const Index& iphase) -> Indices;

/// Get the indices of the species that contains a given element in a multiphase system
/// @param system The multiphase system
/// @param ielement The index of the element in the chemical system
auto indicesSpeciesWithElement(const Multiphase& multiphase, const Index& ielement) -> Indices;

/// Get the index of the phase that contains a given species in a multiphase system
/// @param system The multiphase system
/// @param ispecies The index of the species in the chemical system
auto indexPhaseWithSpecies(const Multiphase& multiphase, const Index& ispecies) -> Index;

/// Get the local index (i.e., the index in its corresponding phase) of a species in a multiphase system
/// @param system The multiphase system
/// @param ispecies The index of the species
auto localIndexSpecies(const Multiphase& multiphase, const Index& ispecies) -> Index;

/// Calculate the standard chemical potentials of the species in a multiphase system
/// @param system The multiphase system
/// @param T The temperature of the system (in units of K)
/// @param P The pressure of the system (in units of Pa)
auto chemicalPotentials(const Multiphase& multiphase, double T, double P) -> Vector;

/// Calculate the molar fractions of the species in a multiphase system
/// @param system The multiphase system
/// @param n The molar amounts of the species (in units of mol)
auto molarFractions(const Multiphase& multiphase, const Vector& n) -> Vector;

/// Calculate the concentrations of the species in a multiphase system
/// @param system The multiphase system
/// @param n The molar amounts of the species (in units of mol)
auto concentrations(const Multiphase& multiphase, const Vector& n) -> Vector;

/// Calculate the activity of a species in a multiphase system
/// @param system The multiphase system
/// @param ispecies The index of the species in the chemical system
/// @param T The temperature of the system (in units of K)
/// @param P The pressure of the system (in units of Pa)
/// @param n The molar amounts of the species (in units of mol)
auto activity(const Multiphase& multiphase, const Index& ispecies, double T, double P, const Vector& n) -> ScalarResult;

/// Calculate the activities of the species in a multiphase system
/// @param system The multiphase system
/// @param T The temperature of the system (in units of K)
/// @param P The pressure of the system (in units of Pa)
/// @param n The molar amounts of the species (in units of mol)
auto activities(const Multiphase& multiphase, double T, double P, const Vector& n) -> VectorResult;

/// Assemble the formula matrix of a multiphase system
///
/// The formula matrix of a multiphase system is defined as a matrix
/// whose entry `(i, j)` is given by the number of atoms of the `i`-th
/// element in the `j`-th species.
///
/// @param system The multiphase system
auto formulaMatrix(const Multiphase& multiphase) -> Matrix;

/// Get the view of a vector that corresponds to the entries of a given phase
/// @param system The multiphase system
/// @param iphase The index of the phase in the chemical system
/// @param vec The vector instance
auto subvector(const Multiphase& multiphase, const Index& iphase, const Vector& vec) -> SubVector;

} /* namespace Reaktor */
