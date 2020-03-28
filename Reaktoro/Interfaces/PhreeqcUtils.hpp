// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// C++ includes
#include <map>
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>

// Phreeqc includes
#include <Reaktoro/Interfaces/PhreeqcLegacy.hpp>

namespace Reaktoro {

namespace PhreeqcUtils {

/// Load the PHREEQC instance with a database file.
/// @param phreeqc The PHREEQC instance
/// @param database The path to the database file
auto load(PHREEQC& phreeqc, std::string database) -> void;

/// Execute a PHREEQC input script file.
/// @param phreeqc The PHREEQC instance
/// @param input The input script either as a file name or as a input string
/// @param output The file name where the result is output
auto execute(PHREEQC& phreeqc, std::string input, std::string output) -> void;

/// Find an element in a Phreeqc instance.
/// @param phreeqc The Phreeqc instance
/// @param name The name of the element
/// @return A pointer to the element if found, nullptr otherwise.
auto findElement(const PHREEQC& phreeqc, std::string name) -> PhreeqcElement*;

/// Find a species in a Phreeqc instance.
/// @param phreeqc The Phreeqc instance
/// @param name The name of the species
/// @return A pointer to the species if found, nullptr otherwise.
auto findSpecies(const PHREEQC& phreeqc, std::string name) -> PhreeqcSpecies*;

/// Find a phase in a Phreeqc instance.
/// @param phreeqc The Phreeqc instance
/// @param name The name of the phase
/// @return A pointer to the phase if found, nullptr otherwise.
auto findPhase(const PHREEQC& phreeqc, std::string name) -> PhreeqcPhase*;

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

/// Return the name of the Phreeqc element according to Reaktoro's convention.
auto name(const PhreeqcElement* element) -> std::string;

/// Return the name of the Phreeqc species according to Reaktoro's convention.
auto name(const PhreeqcSpecies* species) -> std::string;

/// Return the name of the Phreeqc phase according to Reaktoro's convention.
auto name(const PhreeqcPhase* phase) -> std::string;

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

/// Return true if the Phreeqc species instance is an aqueous species.
auto isAqueousSpecies(const PhreeqcSpecies* species) -> bool;

/// Return true if the Phreeqc species instance is an exchange species.
auto isExchangeSpecies(const PhreeqcSpecies* species) -> bool;

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

/// Return the active aqueous species in a Phreeqc instance.
/// @param phreeqc The Phreeqc instance
auto activeAqueousSpecies(const PHREEQC& phreeqc) -> std::vector<PhreeqcSpecies*>;

/// Return the active exchange species in a PHREEQC instance.
auto activeExchangeSpecies(const PHREEQC& phreeqc) -> std::vector<PhreeqcSpecies*>;

/// Return the active aqueous product species in a Phreeqc instance.
/// A product aqueous species is an aqueous species defined in terms of master species.
/// @param phreeqc The Phreeqc instance
auto activeProductSpecies(const PHREEQC& phreeqc) -> std::vector<PhreeqcSpecies*>;

/// Return the active gaseous species in a Phreeqc instance defined in a `GAS_PHASE` block.
/// @param phreeqc The Phreeqc instance
auto activeGaseousSpecies(const PHREEQC& phreeqc) -> std::vector<PhreeqcPhase*>;

/// Return the active phases in a Phreeqc instance defined in a `EQUILIBRIUM_PHASES` block.
/// @param phreeqc The Phreeqc instance
auto activePhasesInEquilibriumPhases(const PHREEQC& phreeqc) -> std::vector<PhreeqcPhase*>;

/// Return the phases in a Phreeqc instance that marked for saturation index calculation.
/// @param phreeqc The Phreeqc instance
auto activePhasesInSaturationList(const PHREEQC& phreeqc) -> std::vector<PhreeqcPhase*>;

/// Return the molar amounts of Phreeqc species (aqueous species)
/// @param phreeqc The Phreeqc instance
/// @param species The container with pointers to Phreeqc  species instances
auto speciesAmounts(const PHREEQC& phreeqc, const std::vector<PhreeqcSpecies*>& species) -> VectorXr;

/// Return the molar amounts of Phreeqc phases (gaseous or mineral species)
/// @param phreeqc The Phreeqc instance
/// @param species The container with pointers to Phreeqc phase instances
auto speciesAmounts(const PHREEQC& phreeqc, const std::vector<PhreeqcPhase*>& phases) -> VectorXr;

/// Return the natural logarithm of the equilibrium constant of a Phreeqc species (aqueous species)
/// @param sspecies A pointer to the Phreeqc species (aqueous species)
/// @param T The temperature (in units of K)
/// @param P The pressure (in units of Pa)
auto lnEquilibriumConstant(const PhreeqcSpecies* species, double T, double P) -> real;

/// Return the natural logarithm of the equilibrium constant of a Phreeqc phase (gaseous or mineral species)
/// @param phase A pointer to the Phreeqc phase (gaseous or mineral species)
/// @param T The temperature (in units of K)
/// @param P The pressure (in units of Pa)
auto lnEquilibriumConstant(const PhreeqcPhase* phase, double T, double P) -> real;

} // namespace PhreeqcUtils
} // namespace Reaktoro
