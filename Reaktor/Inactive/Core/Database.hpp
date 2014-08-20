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
#include <memory>
#include <string>
#include <vector>

namespace Reaktor {

// Reaktor forward declarations
class AqueousSpecies;
class GaseousSpecies;
class MineralSpecies;

/**
 * Provides operations to retrive physical and thermodynamic data of chemical species
 *
 * The Database class is used to retrieve information of chemical species. It is initialized
 * from a `xml` database file, which must satisfy some format conditions. Once it is initialized,
 * one can, for example, retrieve thermodynamic information of an aqueous species that will be
 * used to calculate its standard chemical potential.
 *
 * **Usage**
 *
 * In the example below, a Database instance is initialized and
 * two queries are made to retrieve information of an aqueous and
 * a gaseous species. Note that if a species is not present in the
 * database, then an exception is thrown.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * using namespace Reaktor;
 *
 * // Create a Database instance by parsing a local database file
 * Database database("geodb.xml")
 *
 * // Retrieve information of species H2O(l) and CO2(g)
 * AqueousSpecies aqueousSpecies = database.aqueousSpecies("H2O(l)");
 * GaseousSpecies gaseousSpecies = database.gaseousSpecies("CO2(g)");
 *
 * // Output the data of the species H2O(l) and CO2(g)
 * std::cout << aqueousSpecies << std::endl;
 * std::cout << gaseousSpecies << std::endl;
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * @see AqueousSpecies, GaseousSpecies, MineralSpecies
 * @ingroup Core
 */
class Database
{
public:
    /**
     * Constructs a default Database instance
     */
    Database();

    /**
     * Constructs a Database instance by parsing a `xml` database file
     * @param filename The name of the database file
     */
    explicit Database(const std::string& filename);

    /**
     * Gets all aqueous species in the database
     */
    auto aqueousSpecies() -> std::vector<AqueousSpecies>;

    /**
     * Retrieves information of a aqueous species in the database
     *
     * **Note:** An exception is thrown if the database does not contain the species.
     * @param name The name of the aqueous species
     * @return An AqueousSpecies instance containing the species data
     */
    auto aqueousSpecies(const std::string& name) const -> const AqueousSpecies&;

    /**
     * Gets all gaseous species in the database
     */
    auto gaseousSpecies() -> std::vector<GaseousSpecies>;

    /**
     * Retrieves information of a gaseous species in the database
     *
     * **Note:** An exception is thrown if the database does not contain the species.
     * @param name The name of the gaseous species
     * @return A GaseousSpecies instance containing the species data
     */
    auto gaseousSpecies(const std::string& name) const -> const GaseousSpecies&;

    /**
     * Gets all mineral species in the database
     */
    auto mineralSpecies() -> std::vector<MineralSpecies>;

    /**
     * Retrieves information of a mineral species in the database
     *
     * **Note:** An exception is thrown if the database does not contain the species.
     * @param name The name of the mineral species
     * @return A MineralSpecies instance containing the species data
     */
    auto mineralSpecies(const std::string& name) const -> const MineralSpecies&;

    /**
     * Checks if the database contains a given aqueous species
     * @param species The name of the aqueous species
     * @return `true` if the database contains the species
     */
    auto containsAqueousSpecies(const std::string& species) const -> bool;

    /**
     * Checks if the database contains a given gaseous species
     * @param species The name of the gaseous species
     * @return `true` if the database contains the species
     */
    auto containsGaseousSpecies(const std::string& species) const -> bool;

    /**
     * Checks if the database contains a given mineral species
     * @param species The name of the mineral species
     * @return `true` if the database contains the species
     */
    auto containsMineralSpecies(const std::string& species) const -> bool;

    /**
     * Retrieves a list of aqueous species that contains at least one of the specified elements
     * @param elements The list of chemical elements
     * @return A list of names of aqueous species
     */
    auto aqueousSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<std::string>;

    /**
     * Retrieves a list of gaseous species that contains at least one of the specified elements
     * @param elements The list of chemical elements
     * @return A list of names of gaseous species
     */
    auto gaseousSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<std::string>;

    /**
     * Retrieves a list of mineral species that contains at least one of the specified elements
     * @param elements The list of chemical elements
     * @return A list of names of mineral species
     */
    auto mineralSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<std::string>;

private:
    class Impl;

    std::shared_ptr<Impl> pimpl;
};

} /* namespace Reaktor */
