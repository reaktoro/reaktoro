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

namespace Reaktor {

// Forward declarations
class Species;

/// Get the number of elements in a chemical species
/// @param species The chemical species instance
auto numElements(const Species& species) -> unsigned;

/// Get the index of an element that compose a given chemical species
/// @param species The chemical species instance
/// @param element The name of the element in the chemical species
auto indexElement(const Species& species, const std::string& element) -> Index;

/// Check if a chemical species contains a chemical element
/// @param species The chemical species instance
/// @param element The name of the chemical element
auto containsElement(const Species& species, const std::string& element) -> bool;

/// Get the number of atoms of an element in the species
/// @param species The chemical species
/// @param element The name of the chemical element
auto elementAtoms(const Species& species, const std::string& element) -> double;

/// Calculate the standard chemical potential of a chemical species (in units of J/mol)
/// @param species The chemical species
/// @param T The temperature for the calculation (in units of K)
/// @param P The pressure for the calculation (in units of Pa)
auto chemicalPotential(const Species& species, double T, double P) -> double;

/// Get the names of the chemical species in a container of species
/// @param The container of chemical species
auto names(const std::vector<Species>& species) -> std::vector<std::string>;

/// Get the electrical charges of the chemical species in a container of species
/// @param The container of chemical species
auto charges(const std::vector<Species>& species) -> std::vector<double>;

} // namespace Reaktor
