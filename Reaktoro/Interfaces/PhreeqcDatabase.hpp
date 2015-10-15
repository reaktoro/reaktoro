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
    /// @return A pointer to the Phreeqc element instance if found, nullptr otherwise.
    auto element(std::string name) -> PhreeqcElement*;

    /// Return a pointer to the aqueous species with given name.
    /// @param name The name of the element
    /// @return A pointer to the Phreeqc species instance if found, nullptr otherwise.
    auto species(std::string name) -> PhreeqcSpecies*;

    /// Return a pointer to the gaseous or mineral species with given name.
    /// @param name The name of the phase
    /// @return A pointer to the Phreeqc phase instance if found, nullptr otherwise.
    auto phase(std::string name) -> PhreeqcPhase*;

    /// Return the elements that compose an aqueous species
    /// @param species A pointer to a Phreeqc species
    auto elements(const PhreeqcSpecies* species) -> std::map<PhreeqcElement*, double>;

    /// Return the elements that compose a gaseous or mineral species
    /// @param phase A pointer to a Phreeqc phase
    auto elements(const PhreeqcPhase* phase) -> std::map<PhreeqcElement*, double>;

    /// Return the stoichiometry of an element in a Phreeqc species (aqueous species).
    /// The element name `Z` denotes the charge element.
    /// @param element The name of the element
    /// @param species A pointer to a Phreeqc species (aqueous species)
    auto stoichiometry(std::string element, const PhreeqcSpecies* species) -> double;

    /// Return the stoichiometry of an element a Phreeqc phase (gaseous or mineral species).
    /// The element name `Z` denotes the charge element.
    /// @param element The name of the element
    /// @param phase A pointer to a Phreeqc phase (gaseous or mineral species)
    auto stoichiometry(std::string element, const PhreeqcPhase* phase) -> double;

    /// Return the reaction equation of a Phreeqc species (aqueous species).
    /// The equation is defined by a set of pairs with the names of the species
    /// defining the reaction and their stoichiometry coefficients.
    /// An empty equation is returned in case the given species is a primary species.
    /// @param species A pointer to a Phreeqc species (aqueous species)
    auto reactionEquation(const PhreeqcSpecies* species) -> std::map<std::string, double>;

    /// Return the reaction equation of a Phreeqc phase (gaseous or mineral species).
    /// The equation is defined by a set of pairs with the names of the species
    /// defining the reaction and their stoichiometry coefficients.
    /// @param phase A pointer to a Phreeqc phase (gaseous or mineral species)
    auto reactionEquation(const PhreeqcPhase* phase) -> std::map<std::string, double>;

    /// Return the aqueous species in the Phreeqc database.
    auto aqueousSpecies() -> std::vector<PhreeqcSpecies*>;

    /// Return the gaseous species in the Phreeqc database.
    auto gaseousSpecies() -> std::vector<PhreeqcPhase*>;

    /// Return the mineral species in the Phreeqc database.
    auto mineralSpecies() -> std::vector<PhreeqcPhase*>;

	/// Return the master aqueous species in the Phreeqc database.
    /// The master species are those species that serve as building blocks for secondary species.
    auto masterSpecies() -> std::vector<PhreeqcSpecies*>;

    /// Return the secondary aqueous species in the Phreeqc database.
    /// The secondary aqueous species are those species constructed from master species.
    auto secondarySpecies() -> std::vector<PhreeqcSpecies*>;

    /// Return the ln equilibrium constant of a Phreeqc species (aqueous species)
    /// @param species A pointer to the Phreeqc species (aqueous species)
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    auto lnEquilibriumConstant(const PhreeqcSpecies* species, double T, double P) -> double;

    /// Return the ln equilibrium constant of a Phreeqc phase (gaseous or mineral species)
    /// @param phase A pointer to the Phreeqc phase (gaseous or mineral species)
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    auto lnEquilibriumConstant(const PhreeqcPhase* phase, double T, double P) -> double;
};

} // namespace Reaktoro
