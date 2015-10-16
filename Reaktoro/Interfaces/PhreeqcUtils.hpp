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
#include <map>
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>

// Phreeqc includes
#include <Reaktoro/Interfaces/PHREEQC.hpp>

namespace Reaktoro {
namespace PhreeqcUtils {

/// Return the elements that compose a Phreeqc species (aqueous species)
/// @param sspecies A pointer to a Phreeqc species
auto elements(const PhreeqcSpecies* species) -> std::map<PhreeqcElement*, double>;

/// Return the elements that compose a Phreeqc phase (gaseous or mineral species)
/// @param phase A pointer to a Phreeqc phase
auto elements(const PhreeqcPhase* phase) -> std::map<PhreeqcElement*, double>;

/// Return the stoichiometry of an element in a Phreeqc species (aqueous species).
/// The element name `Z` denotes the charge element.
/// @param element The name of the element
/// @param sspecies A pointer to a Phreeqc species (aqueous species)
auto stoichiometry(std::string element, const PhreeqcSpecies* species) -> double;

/// Return the stoichiometry of an element in a Phreeqc phase (gaseous or mineral species).
/// The element name `Z` denotes the charge element.
/// @param element The name of the element
/// @param phase A pointer to a Phreeqc phase (gaseous or mineral species)
auto stoichiometry(std::string element, const PhreeqcPhase* phase) -> double;

/// Return the reaction equation of a Phreeqc species (aqueous species).
/// The equation is defined by a map of the names of the species
/// defining the reaction and their stoichiometry coefficients.
/// An empty equation is returned in case the given species is a primary species.
/// @param sspecies A pointer to a Phreeqc species (aqueous species)
auto reactionEquation(const PhreeqcSpecies* species) -> std::map<std::string, double>;

/// Return the reaction equation of a Phreeqc phase (gaseous or mineral species).
/// The equation is defined by a map of the names of the species
/// defining the reaction and their stoichiometry coefficients.
/// @param phase A pointer to a Phreeqc phase (gaseous or mineral species)
auto reactionEquation(const PhreeqcPhase* phase) -> std::map<std::string, double>;

/// Return true if the Phreeqc phase instance is a gaseous species.
auto isGaseousSpecies(const PhreeqcPhase* phase) -> bool;

/// Return true if the Phreeqc phase instance is a mineral species.
auto isMineralSpecies(const PhreeqcPhase* phase) -> bool;

/// Return the index of a Phreeqc species (aqueous species) in a set of species.
/// @param name The name of the Phreeqc species
/// @param species The container with pointers to Phreeqc  species instances
auto index(std::string name, const std::vector<PhreeqcSpecies*>& species) -> unsigned;

/// Return the index of a Phreeqc phase (gaseous or mineral species) in a set of phases.
/// @param name The name of the Phreeqc phase
/// @param phases The container with pointers to Phreeqc phase instances
auto index(std::string name, const std::vector<PhreeqcPhase*>& phases) -> unsigned;

/// Return the molar amounts of Phreeqc species (aqueous species)
/// @param species The container with pointers to Phreeqc  species instances
auto speciesAmounts(const std::vector<PhreeqcSpecies*>& species) -> Vector;

/// Return the molar amounts of Phreeqc phases (gaseous or mineral species)
/// @param species The container with pointers to Phreeqc phase instances
auto speciesAmounts(const std::vector<PhreeqcPhase*>& phases) -> Vector;

/// Return the natural logarithm of the equilibrium constant of a Phreeqc species (aqueous species)
/// @param sspecies A pointer to the Phreeqc species (aqueous species)
/// @param T The temperature (in units of K)
/// @param P The pressure (in units of Pa)
auto lnEquilibriumConstant(const PhreeqcSpecies* species, double T, double P) -> double;

/// Return the natural logarithm of the equilibrium constant of a Phreeqc phase (gaseous or mineral species)
/// @param phase A pointer to the Phreeqc phase (gaseous or mineral species)
/// @param T The temperature (in units of K)
/// @param P The pressure (in units of Pa)
auto lnEquilibriumConstant(const PhreeqcPhase* phase, double T, double P) -> double;

} // namespace PhreeqcUtils
} // namespace Reaktoro
