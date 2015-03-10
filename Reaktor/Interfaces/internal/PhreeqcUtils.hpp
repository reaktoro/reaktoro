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
#include <map>
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Common/Matrix.hpp>

// Forward declarations of Phreeqc types
class Phreeqc;
class element;
class species;
class phase;

namespace Reaktor {

/// Load a database in Phreeqc
/// @param phreeqc The Phreeqc instance
/// @param filename The name of the database file
auto loadDatabase(Phreeqc& phreeqc, std::string filename) -> void;

/// Load and execute a Phreeqc script
/// @param phreeqc The Phreeqc instance
/// @param filename The name of the script file
auto loadScript(Phreeqc& phreeqc, std::string filename) -> void;

/// Find an element in a Phreeqc instance
/// @param phreeqc The Phreeqc instance
/// @param name The name of element
/// @return A pointer to the element if found, nullptr otherwise.
auto findElement(Phreeqc& phreeqc, std::string name) -> element*;

/// Find an aqueous species in a Phreeqc instance
/// @param phreeqc The Phreeqc instance
/// @param name The name of the aqueous species
/// @return A pointer to the species if found, nullptr otherwise.
auto findSpecies(Phreeqc& phreeqc, std::string name) -> species*;

/// Find a phase (a gas or mineral species) in a Phreeqc instance
/// @param phreeqc The Phreeqc instance
/// @param name The name of the phase
/// @return A pointer to the phase if found, nullptr otherwise.
auto findPhase(Phreeqc& phreeqc, std::string name) -> phase*;

/// Get the elements that compose a Phreeqc aqueous species
/// @param s A pointer to a Phreeqc species
/// @return A map of element pointers and their coefficients.
auto getElementsInSpecies(const species* s) -> std::map<element*, double>;

/// Get the elements that compose a Phreeqc phase (a gas or mineral)
/// @param s A pointer to a Phreeqc species
/// @return A map of element pointers and their coefficients.
auto getElementsInPhase(const phase* p) -> std::map<element*, double>;

/// Get the reaction equation of a Phreeqc species (aqueous species).
/// The equation is defined by a map of the names of the species
/// defining the reaction and their stoichiometry coefficients.
/// An empty equation is returned in case the given species is a primary species.
/// @param s A pointer to a Phreeqc species (aqueous species)
auto getReactionEquation(const species* s) -> std::map<std::string, double>;

/// Get the reaction equation of a Phreeqc phase (gas or mineral).
/// The equation is defined by a map of the names of the species
/// defining the reaction and their stoichiometry coefficients.
/// @param p A pointer to a Phreeqc phase (gas or mineral)
auto getReactionEquation(const phase* p) -> std::map<std::string, double>;

/// Get the number of element atoms in a Phreeqc species (aqueous species).
/// @param element The name of the element
/// @param s A pointer to a Phreeqc species (aqueous species)
auto numElementAtomsInSpecies(std::string element, const species* s) -> double;

/// Get the number of element atoms in a Phreeqc phase (gas or mineral).
/// @param element The name of the element
/// @param s A pointer to a Phreeqc species (gas or mineral)
auto numElementAtomsInPhase(std::string element, const phase* p) -> double;

/// Collect the active aqueous species in a Phreeqc instance.
/// @param phreeqc The Phreeqc instance
auto collectAqueousSpecies(Phreeqc& phreeqc) -> std::vector<species*>;

/// Collect the active secondary aqueous species in a Phreeqc instance.
/// @param phreeqc The Phreeqc instance
auto collectSecondarySpecies(Phreeqc& phreeqc) -> std::vector<species*>;

/// Collect the active gaseous species in a Phreeqc instance.
/// @param phreeqc The Phreeqc instance
auto collectGaseousSpecies(Phreeqc& phreeqc) -> std::vector<phase*>;

/// Collect the active mineral species in a Phreeqc instance.
/// @param phreeqc The Phreeqc instance
auto collectMineralSpecies(Phreeqc& phreeqc) -> std::vector<phase*>;

/// Get the index of a Phreeqc species (aqueous species) in a set of species.
/// @param name The name of the Phreeqc species
/// @param pointers The container of Phreeqc species pointers
auto index(std::string name, const std::vector<species*>& pointers) -> unsigned;

/// Get the index of a Phreeqc phase (gas or mineral) in a set of phases.
/// @param name The name of the Phreeqc phase
/// @param pointers The container of Phreeqc phase pointers
auto index(std::string name, const std::vector<phase*>& pointers) -> unsigned;

/// Return the molar amounts of Phreeqc species (aqueous species)
/// @param pointers The container of Phreeqc species
auto speciesAmountsInSpecies(const std::vector<species*>& pointers) -> Vector;

/// Return the molar amounts of Phreeqc phases (gas or mineral)
/// @param pointers The container of Phreeqc phases
auto speciesAmountsInPhases(const std::vector<phase*>& pointers) -> Vector;

/// Return the natural logarithm of the equilibrium constant of a Phreeqc species (aqueous species)
/// @param s A pointer to the Phreeqc species (aqueous species)
/// @param T The temperature (in units of K)
/// @param P The pressure (in units of Pa)
auto lnEquilibriumConstant(const species* s, double T, double P) -> double;

/// Return the natural logarithm of the equilibrium constant of a Phreeqc phase (gas or mineral)
/// @param p A pointer to the Phreeqc phase (gas or mineral)
/// @param T The temperature (in units of K)
/// @param P The pressure (in units of Pa)
auto lnEquilibriumConstant(const phase* p, double T, double P) -> double;

} // namespace Reaktor
