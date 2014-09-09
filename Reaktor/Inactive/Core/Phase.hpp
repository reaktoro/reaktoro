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
#include <memory>
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/Matrix.hpp>
#include <Reaktor/Common/Vector.hpp>
#include <Reaktor/Core/Types.hpp>

namespace Reaktor {

// Forward declarations
class Species;

/**
 * Provides a computational representation of a phase
 *
 * The Phase class can be seen as a collection of Species instances with
 * functionalities to compute their standard chemical potentials and activities.
 *
 * Instances of class Phase are the building blocks of a ChemicalSystem instance,
 * which computationaly represents multiphase system.
 *
 * A Phase instance can be produced from an AqueousPhase, GaseousPhase, or
 * MineralPhase instance. These classes should be used together with the Database
 * class to define a phase -- specifying their species, their standard chemical
 * potential models, and their activity models.
 *
 * @see AqueousPhase, GaseousPhase, MineralPhase, Database, ChemicalSystem
 * @ingroup Core
 */
class Phase
{
public:
    /**
     * Construct a default Phase instance
     */
    Phase();

    /**
     * Constructs a copy of a Phase instance
     */
    Phase(const Phase& other);

    /**
     * Destroys this phase instance
     */
    virtual ~Phase();

    /**
     * Assigns a Phase instance to this instance
     */
    auto operator=(Phase other) -> Phase&;

    /**
     * Sets the name of the phase
     * @param name The name of the phase
     * @return A reference to this Phase instance
     */
    auto setName(const std::string& name) -> Phase&;

    /**
     * Sets the chemical species that compose the phase
     * @param species The species that compose the phase
     * @return A reference to this Phase instance
     * @see Species
     */
    auto setSpecies(const std::vector<Species>& species) -> Phase&;

    /**
     * Sets the concentration function of the phase
     * @param concentration The concentration function of the phase
     * @return A reference to this Phase instance
     * @see ConcentrationFn
     */
    auto setConcentrationFn(const ConcentrationFn& concentration) -> Phase&;

    /**
     * Sets the activity function of the phase
     * @param activity The activity function
     * @return A reference to this Phase instance
     * @see Activity
     */
    auto setActivityFn(const ActivityFn& activity) -> Phase&;

    /**
     * Gets the name of the phase
     */
    auto name() const -> const std::string&;

    /**
     * Gets the number of species in the phase
     */
    auto numSpecies() const -> unsigned;

    /**
     * Gets the index of a species in the phase
     * @param species The name of the species
     * @return The index of the species if found. The number of species otherwise.
     */
    auto idxSpecies(const std::string& species) const -> Index;

    /**
     * Gets a species in the phase
     * @param idxSpecies The index of the species
     * @return A const reference to the species
     */
    auto species(const Index& idxSpecies) const -> const Species&;

    /**
     * Gets a species in the phase
     * @param idxSpecies The index of the species
     * @return A reference to the species
     */
    auto species(const Index& idxSpecies) -> Species&;

    /**
     * Gets a species in the phase
     * @param species The name of the species
     * @return A const reference to the species
     */
    auto species(const std::string& species) const -> const Species&;

    /**
     * Gets a species in the phase
     * @param species The name of the species
     * @return A reference to the species
     */
    auto species(const std::string& species) -> Species&;

    /**
     * Gets the species that compose the phase
     * @return A const reference to the species
     */
    auto species() const -> const std::vector<Species>&;

    /**
     * Gets the species that compose the phase
     * @return A reference to the species
     */
    auto species() -> std::vector<Species>&;

    /**
     * Gets the names of the species that compose the phase
     */
    auto speciesNames() const -> std::vector<std::string>;

    /**
     * Calculates the standard chemical potentials of the species (in units of J/mol)
     * @param T The temperature of the chemical system (in units of K)
     * @param P The pressure of the chemical system (in units of Pa)
     * @return A vector of standard chemical potentials (in units of J/mol)
     */
    auto chemicalPotentials(double T, double P) const -> Vector;

    /**
     * Calculates the molar fractions of the species
     * @param n The molar abundance of the species (in units of mol)
     * @return The molar fractions of the species
     */
    auto molarFractions(const Vector& n) const -> Vector;

    /**
     * Calculates the concentrations of the species
     * @param n The molar abundance of the species (in units of mol)
     * @return The concentrations of the species
     */
    auto concentrations(const Vector& n) const -> Vector;

    /**
     * Calculates the activities and their molar derivatives of the species
     * @param T The temperature of the chemical system (in units of K)
     * @param P The pressure of the chemical system (in units of Pa)
     * @param n The molar abundance of the species (in units of mol)
     * @return A activities and their molar derivatives of the species
     */
    auto activities(double T, double P, const Vector& n) const -> VectorResult;

    /**
     * Checks if this Phase instance is equal another
     */
    auto operator==(const Phase& phase) const -> bool;

private:
    class Impl;

    std::unique_ptr<Impl> pimpl;
};

/**
 * Outputs a Phase instance
 */
auto operator<<(std::ostream& out, const Phase& phase) -> std::ostream&;

} // namespace Reaktor
