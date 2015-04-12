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
#include <string>
#include <vector>
#include <map>

namespace Reaktoro {

/// Return a vector of all known 116 chemical elements.
auto elements() -> std::vector<std::string>;

/// Determine. the elemental composition of a chemical compound.
/// @param compound The formula of the chemical compound
/// @return The elemental composition of the chemical compound
auto elements(std::string formula) -> std::map<std::string, double>;

/// Return the atomic mass of a chemical element (in units of kg/mol).
/// @param element The symbol of the chemical element
/// @return The atomic mass of the chemical element
auto atomicMass(std::string element) -> double;

/// Calculate the molar mass of a chemical compound (in units of kg/mol).
/// @param compound The elements and their stoichiometries that compose the chemical compound
/// @return The molar mass of the chemical compound
auto molarMass(const std::map<std::string, double>& compound) -> double;

/// Calculate the molar mass of a chemical compound (in units of kg/mol).
/// @param compound The formula of the chemical compound
/// @return The molar mass of the chemical compound
auto molarMass(std::string formula) -> double;

/// Extract the electrical charge of a chemical formula.
/// The chemical formula of a species can be represented by
/// string with the following formats:
/// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// std::string formula1 = "Na+";   // species: Na+
/// std::string formula2 = "Ca+2";  // species: Ca++
/// std::string formula3 = "OH-";   // species: OH-
/// std::string formula4 = "CO3-2"; // species: CO3--
/// std::string formula5 = "H+1";   // species: H+
/// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// The number 1 is optional when the species has one
/// negative or positive electrical charge.
/// @param formula The chemical formula of a species
/// @return The electrical charge of the species
auto charge(const std::string& formula) -> double;

} // namespace Reaktoro
