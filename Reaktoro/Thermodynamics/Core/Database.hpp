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
#include <memory>
#include <string>
#include <vector>

namespace Reaktoro {

// Forward declarations
class Element;
class AqueousSpecies;
class GaseousSpecies;
class MineralSpecies;

/// Provides operations to retrive physical and thermodynamic data of chemical species.
///
/// The Database class is used to retrieve information of chemical species. It is initialized
/// from a `xml` database file, which must satisfy some format conditions. Once it is initialized,
/// one can, for example, retrieve thermodynamic information of an aqueous species that will be
/// used to calculate its standard chemical potential.
///
//////*Usage**
///
/// In the example below, a Database instance is initialized and
/// two queries are made to retrieve information of an aqueous and
/// a gaseous species. Note that if a species is not present in the
/// database, then an exception is thrown.
/// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// using namespace Reaktoro;
///
/// // Create a Database instance by parsing a local database file
/// Database database("geodb.xml")
///
/// // Retrieve information of species H2O(l) and CO2(g)
/// AqueousSpecies aqueousSpecies = database.aqueousSpecies("H2O(l)");
/// GaseousSpecies gaseousSpecies = database.gaseousSpecies("CO2(g)");
///
/// // Output the data of the species H2O(l) and CO2(g)
/// std::cout << aqueousSpecies << std::endl;
/// std::cout << gaseousSpecies << std::endl;
/// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
///
/// @see AqueousSpecies, GaseousSpecies, MineralSpecies
/// @ingroup Core
class Database
{
public:
    /// Construct a default Database instance
    Database();

    /// Construct a Database instance by parsing a `xml` database file.
    /// If `filename` does not point to a valid database file or the
    /// database file is not found, then a default built-in database
    /// with the same name will be tried. If no default built-in database
    /// exist with given name, an exception will be thrown.
    /// @param filename The name of the database file
    explicit Database(std::string filename);

    /// Add an Element instance in the database.
    auto addElement(const Element& element) -> void;

    /// Add an AqueousSpecies instance in the database.
    auto addAqueousSpecies(const AqueousSpecies& species) -> void;

    /// Add a GaseousSpecies instance in the database.
    auto addGaseousSpecies(const GaseousSpecies& species) -> void;

    /// Add a MineralSpecies instance in the database.
    auto addMineralSpecies(const MineralSpecies& species) -> void;

    /// Return all elements in the database
    auto elements() -> std::vector<Element>;

    /// Return all aqueous species in the database
    auto aqueousSpecies() -> std::vector<AqueousSpecies>;

    /// Return an aqueous species in the database.
    /// **Note:** An exception is thrown if the database does not contain the species.
    /// @param name The name of the aqueous species
    auto aqueousSpecies(std::string name) const -> const AqueousSpecies&;

    /// Return all gaseous species in the database
    auto gaseousSpecies() -> std::vector<GaseousSpecies>;

    /// Return a gaseous species in the database.
    /// **Note:** An exception is thrown if the database does not contain the species.
    /// @param name The name of the gaseous species
    auto gaseousSpecies(std::string name) const -> const GaseousSpecies&;

    /// Return all mineral species in the database
    auto mineralSpecies() -> std::vector<MineralSpecies>;

    /// Return a mineral species in the database.
    /// **Note:** An exception is thrown if the database does not contain the species.
    /// @param name The name of the mineral species
    auto mineralSpecies(std::string name) const -> const MineralSpecies&;

    /// Check if the database contains a given aqueous species
    /// @param species The name of the aqueous species
    auto containsAqueousSpecies(std::string species) const -> bool;

    /// Check if the database contains a given gaseous species
    /// @param species The name of the gaseous species
    auto containsGaseousSpecies(std::string species) const -> bool;

    /// Check if the database contains a given mineral species
    /// @param species The name of the mineral species
    auto containsMineralSpecies(std::string species) const -> bool;

    /// Return the aqueous species that contains at least one of the specified elements.
    auto aqueousSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<AqueousSpecies>;

    /// Return the gaseous species that contains at least one of the specified elements.
    auto gaseousSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<GaseousSpecies>;

    /// Return the mineral species that contains at least one of the specified elements.
    auto mineralSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<MineralSpecies>;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
