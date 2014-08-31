/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

// C++ includes
#include <iostream>
#include <string>
#include <vector>
#include <memory>

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/Matrix.hpp>
#include <Reaktor/Common/PartialVector.hpp>
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

// Forward declarations
class Phase;
class Species;

/**
 * Provides a computational representation of a multiphase chemical system
 *
 * The ChemicalSystem class is used to computationaly model a multiphase chemical system.
 * It is instantiated from a list of Phase instances defining the phases of the system.
 * Once it is built, it cannot be modified (i.e., it is an immutable object).
 *
 * **Usage**
 *
 * A convenient way of building a ChemicalSystem instance is using the ChemicalEditor
 * class. Below we demonstrate how a multiphase system with an aqueous phase, a gaseous
 * phase and three mineral phases can be defined.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * using namespace Reaktor;
 *
 * // Create Database instance from a local database file
 * Database database("geodb.xml");
 *
 * // Create the ChemicalEditor instance to setup the chemical system
 * ChemicalEditor editor(database);
 *
 * // Define an aqueous phase and set the Drummond (1981) activity model for species CO2(aq)
 * editor.addAqueousPhase("H2O(l), H+, OH-, HCO3-, CO2(aq), Ca++, Mg++, CaCl+, MgCl+, CaCl2(aq)")
 *     .setActivityModelDrummondCO2();
 *
 * // Define a gaseous phase and set the Spycher and Pruess (2003) activity model for species H2O(g) and CO2(g)
 * editor.addGaseousPhase("H2O(g), CO2(g)")
 *     .setActivityModelSpycherPruessH2OCO2();
 *
 * // Define some pure mineral phases, with default ideal activity model for the species
 * editor.addMineralPhase("Calcite");
 * editor.addMineralPhase("Dolomite");
 * editor.addMineralPhase("Magnesite");
 *
 * // Create the ChemicalSystem instance
 * ChemicalSystem system = editor;
 *
 * // Output the chemical system
 * std::cout << system << std::endl;
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * @see Phase, Species, ChemicalEditor, Partitioning
 * @ingroup Core
 */
class ChemicalSystem
{
public:
    /**
     * Constructs a ChemicalSystem instance
     */
    ChemicalSystem();

	/**
	 * Constructs a ChemicalSystem instance with a given list of phases
	 * @param phases The list of phases that define the chemical system
	 * @see Phase
	 */
	explicit ChemicalSystem(const std::vector<Phase>& phases);

	/**
	 * Gets the number of phases in the system
	 */
	auto numPhases() const -> unsigned;

	/**
	 * Gets the number of species in the system
	 */
	auto numSpecies() const -> unsigned;

	/**
	 * Gets the number of elements in the system
	 */
	auto numElements() const -> unsigned;

	/**
	 * Gets the phase with given index
	 * @param idx_phase The index of the phase
	 */
	auto phase(const Index& idx_phase) const -> const Phase&;

	/**
	 * Gets the phases in the system
	 */
	auto phases() const -> const std::vector<Phase>&;

	/**
	 * Gets the names of the phases in the system
	 */
	auto phasesNames() const -> std::vector<std::string>;

	/**
     * Gets the chemical species with given index
     * @param idxSpecies The index of the species
     */
    auto species(const Index& idxSpecies) const -> const Species&;

    /**
     * Gets a chemical species with given name
     * @param name The name of the chemical species
     */
    auto species(const std::string& name) const -> const Species&;

	/**
	 * Gets the chemical species in the system
	 */
    auto species() const -> const std::vector<Species>&;

    /**
     * Gets the names of the chemical species in the system
     */
    auto speciesNames() const -> std::vector<std::string>;

    /**
     * Gets the electrical charges of the species
     */
    auto speciesCharges() const -> Vector;

    /**
     * Gets the name of a chemical element in the system
     * @param idx_element The index of the chemical element
     */
    auto element(const Index& idx_element) const -> const std::string&;

	/**
	 * Gets the chemical elements in the system
	 */
	auto elements() const -> const std::vector<std::string>&;

	/**
	 * Gets the formula matrix of the chemical species
	 */
	auto formulaMatrix() const -> const Matrix&;

	/**
	 * Gets the indices of the species in a phase
	 * @param idx_phase The index of the phase
	 * @return The indices of the species in phase with index idxPhase
	 */
    auto idxSpeciesInPhase(const Index& idx_phase) const -> const Indices&;

    /**
     * Gets the indices of the elements in a species
     * @param idx_species The index of the species
     * @return The indices of the elements in the species with index idxSpecies
     */
    auto idxElementsInSpecies(const Index& idx_species) const -> const Indices&;

	/**
	 * Gets the mapping from species to phases
	 *
	 * The mapping *species-to-phases* can be used to determine
	 * the index of the phase where a species belongs to.
	 *
	 * **Usage**
	 *
	 * Below we show how to output every species in the system together with
	 * the name of the phase that contains it.
	 *
	 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	 * // Let system denote a ChemicalSystem instance
	 * const Indices& map = system.mapSpeciesToPhase();
	 *
	 * // Iterate over all species in the system
	 * for(Index i = 0; i < system.numSpecies(); ++i)
	 * {
	 *     const Species& species = system.species(i);
	 *     const Phase& phase = system.phase(map[i]);
	 *
	 *     std::cout << species.name() << " belongs to phase " << phase.name() << std::endl;
	 * }
	 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	 *
	 * @return The indices of the phases of each species
	 */
	auto mapSpeciesToPhase() const -> const Indices&;

	/**
	 * Gets the mapping from phases to species
	 *
	 * The mapping *phases-to-species* can be used to determine
     * the indices of the species that belongs to each phase.
     *
     * **Usage**
     *
     * Below we show how to output every phase in the system together with the
     * names of the species that belongs to it.
     *
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * // Let system denote a ChemicalSystem instance
     * const std::vector<Indices>& map = system.mapPhaseToSpecies();
     *
     * // Iterate over all phases in the system
     * for(Index i = 0; i < system.numPhases(); ++i)
     * {
     *     std::cout << "Phase " << system.phase(i).name() << " contains: ";
     *
     *     // Iterate over all species in the current phase
     *     for(Index j : map[i])
     *         std::cout << system.species(j).name() << " ";
     *
     *     std::cout << std::endl;
     * }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     * @return The indices of the species of each phase
	 */
	auto mapPhaseToSpecies() const -> const std::vector<Indices>&;

	/**
	 * Gets the mapping from species to elements
	 *
	 * The mapping *species-to-elements* can be used to determine
     * the indices of the elements that belongs to each species.
     *
	 * @return The indices of the elements of each species
	 */
	auto mapSpeciesToElements() const -> const std::vector<Indices>&;

	/**
	 * Gets the mapping from elements to species
	 *
	 * The mapping *elements-to-species* can be used to determine
     * the indices of the species that contains a given element.
     *
     * @return The indices of the species that contains each element
	 */
	auto mapElementToSpecies() const -> const std::vector<Indices>&;

	/**
     * Determines the index of an element in the system
     * @param element The name of the element
     * @return The index of the element if it exist. The number of elements otherwise.
     */
    auto idxElement(const std::string& element) const -> Index;

    /**
     * Determines the indices of the elements in the system
     * @param elements The list of element names
     * @return The index of each provided element.
     * If an element does not exist, then its index equals the number of elements.
     */
    auto idxElements(const std::vector<std::string>& elements) const -> Indices;

    /**
     * Determines the index of a species in the system
     * @param species The name of the species
     * @return The index of the species if it exist. The number of species otherwise.
     */
	auto idxSpecies(const std::string& species) const -> Index;

	/**
     * Determines the indices of the species in the system
     * @param species The list of species names
     * @return The index of each provided species.
     * If a species does not exist, then its index equals the number of species.
     */
	auto idxSpecies(const std::vector<std::string>& species) const -> Indices;

	/**
     * Determines the index of a phase in the system
     * @param phase The name of the phase
     * @return The index of the phase if it exist. The number of phases otherwise.
     */
    auto idxPhase(const std::string& phase) const -> Index;

    /**
     * Determines the indices of the phases in the system
     * @param phases The list of phase names
     * @return The index of each provided phase.
     * If an phase does not exist, then its index equals the number of phases.
     */
    auto idxPhases(const std::vector<std::string>& phases) const -> Indices;

	/**
     * Determines the index of the phase that contains a species
     * @param ispecies The index of the species
     * @return The index of the phase containing the species.
     */
    auto idxPhaseWithSpecies(Index ispecies) const -> Index;

    /**
     * Determines the index of the phase that contains a species
     * @param species The name of the species
     * @return The index of the phase containing the species if it exist. The number of phases otherwise.
     */
    auto idxPhaseWithSpecies(const std::string& species) const -> Index;

    /**
     * Determines if a species is contained in the system
     * @param species The name of the species
     * @return `true` if the system contains the species
     */
    auto containsSpecies(const std::string& species) const -> bool;

    /**
     * Determines if a element is contained in the system
     * @param element The name of the element
     * @return `true` if the system contains the element
     */
    auto containsElement(const std::string& element) const -> bool;

    /**
     * Determines if a phase is contained in the system
     * @param phase The name of the phase
     * @return `true` if the system contains the phase
     */
    auto containsPhase(const std::string& phase) const -> bool;

    /**
     * Determines the index of an element in the system and raises an exception if it is not contained
     * @param element The name of the element
     */
    auto idxElementWithError(const std::string& element) const -> Index;

    /**
     * Determines the index of a species in the system and raises an exception if it is not contained
     * @param species The name of the species
     */
    auto idxSpeciesWithError(const std::string& species) const -> Index;

    /**
     * Determines the index of a phase in the system and raises an exception if it is not contained
     * @param phase The name of the phase
     */
    auto idxPhaseWithError(const std::string& phase) const -> Index;

	/**
	 * Calculates the standard chemical potentials of the species
	 * @param T The temperature of the system (in units of K)
	 * @param P The pressure of the system (in units of Pa)
	 */
    auto chemicalPotentials(double T, double P) const -> Vector;

    /**
     * Calculates the molar fractions of the species
     * @param n The molar abundance of the species (in units of mol)
     */
    auto molarFractions(const Vector& n) const -> Vector;

    /**
     * Calculates the concentrations of the species
     * @param n The molar abundance of the species (in units of mol)
     */
    auto concentrations(const Vector& n) const -> Vector;

	/**
	 * Calculates the activities of the species
	 * @param T The temperature of the system (in units of K)
     * @param P The pressure of the system (in units of Pa)
	 * @param n The molar abundance of the species (in units of mol)
	 */
    auto activities(double T, double P, const Vector& n) const -> PartialVector;

private:
    class Impl;

    std::shared_ptr<Impl> pimpl;
};

/**
 * Outputs a ChemicalSystem instance
 */
auto operator<<(std::ostream& out, const ChemicalSystem& system) -> std::ostream&;

} // namespace Reaktor
