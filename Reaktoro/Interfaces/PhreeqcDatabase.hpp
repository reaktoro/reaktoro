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

// Phreeqc includes
#include <Reaktoro/Interfaces/PHREEQC.hpp>

namespace Reaktoro {

class PhreeqcDatabase
{
public:
    /// Construct a default PhreeqcDatabase instance
    PhreeqcDatabase();

    /// Construct a custom PhreeqcDatabase instance
    /// @param filename The path to the Phreeqc database file
    explicit PhreeqcDatabase(std::string filename);

    /// Return a pointer to the element with given name.
    /// @param name The name of the element
    /// @return A pointer to the element instance if found, nullptr otherwise.
    auto element(std::string name) -> element*;

    /// Return a pointer to the aqueous species with given name.
    /// @param name The name of the element
    /// @return A pointer to the species instance if found, nullptr otherwise.
    auto species(std::string name) -> species*;

    /// Return a pointer to the gaseous or mineral species with given name.
    /// @param name The name of the phase
    /// @return A pointer to the phase instance if found, nullptr otherwise.
    auto phase(std::string name) -> phase*;

    /// Return the elements that compose a Phreeqc aqueous species
    /// @param s A pointer to a Phreeqc species
    auto getElements(const species* s) -> std::set<element*>;

    /// Return the elements that compose a Phreeqc phase (a gas or mineral)
    /// @param p A pointer to a Phreeqc phase
    auto getElements(const phase* p) -> std::set<element*>;

    /// Return the stoichiometry of an element in a Phreeqc species (aqueous species).
    /// Note that if an element name with valence is provided, e.g. `N(-3)`,
    /// then the element name without valence, e.g. `N`, is checked instead.
    /// The element name `Z` denotes the charge element.
    /// @param element The name of the element
    /// @param s A pointer to a Phreeqc species (aqueous species)
    auto getElementStoichiometry(std::string element, const species* s) -> double;

    /// Return the stoichiometry of an element a Phreeqc phase (gas or mineral).
    /// Note that if an element name with valence is provided, e.g. `N(-3)`,
    /// then the element name without valence, e.g. `N`, is checked instead.
    /// The element name `Z` denotes the charge element.
    /// @param element The name of the element
    /// @param p A pointer to a Phreeqc phase (gas or mineral)
    auto getElementStoichiometry(std::string element, const phase* p) -> double;

    /// Return the reaction equation of a Phreeqc species (aqueous species).
    /// The equation is defined by a map of the names of the species
    /// defining the reaction and their stoichiometry coefficients.
    /// An empty equation is returned in case the given species is a primary species.
    /// @param s A pointer to a Phreeqc species (aqueous species)
    auto getReactionEquation(const species* s) -> std::map<std::string, double>;

    /// Return the reaction equation of a Phreeqc phase (gas or mineral).
    /// The equation is defined by a map of the names of the species
    /// defining the reaction and their stoichiometry coefficients.
    /// @param p A pointer to a Phreeqc phase (gas or mineral)
    auto getReactionEquation(const phase* p) -> std::map<std::string, double>;

    /// Collect the active aqueous species in a Phreeqc instance.
    /// @param phreeqc The Phreeqc instance
    auto collectAqueousSpecies(PHREEQC& phreeqc) -> std::vector<species*>;

    /// Collect the active secondary aqueous species in a Phreeqc instance.
    /// @param phreeqc The Phreeqc instance
    auto collectSecondarySpecies(PHREEQC& phreeqc) -> std::vector<species*>;

    /// Collect the active gaseous species in a Phreeqc instance.
    /// The collected gaseous species are those defined in a `GAS_PHASE` block.
    /// @param phreeqc The Phreeqc instance
    auto collectGaseousSpecies(PHREEQC& phreeqc) -> std::vector<phase*>;

    /// Collect the gaseous species in the speciation list of a Phreeqc instance.
    /// @param phreeqc The Phreeqc instance
    auto collectGaseousSpeciesInSpeciationList(PHREEQC& phreeqc) -> std::vector<phase*>;

    /// Collect the equilibrium mineral species in a Phreeqc instance.
    /// The collected mineral species are those defined in a `EQUILIBRIUM_PHASES` block.
    /// @param phreeqc The Phreeqc instance
    auto collectMineralSpecies(PHREEQC& phreeqc) -> std::vector<phase*>;

    /// Collect the mineral species in the speciation list of a Phreeqc instance.
    /// @param phreeqc The Phreeqc instance
    auto collectGaseousSpeciesInSpeciationList(PHREEQC& phreeqc) -> std::vector<phase*>;

    /// Return the index of a Phreeqc species (aqueous species) in a set of species.
    /// @param name The name of the Phreeqc species
    /// @param pointers The container of Phreeqc species pointers
    auto index(std::string name, const std::vector<species*>& pointers) -> unsigned;

    /// Return the index of a Phreeqc phase (gas or mineral) in a set of phases.
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


};

} // namespace Reaktoro
