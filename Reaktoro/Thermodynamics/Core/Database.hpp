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


// Forward declarations for ThermoFun
namespace ThermoFun {

class Database;

} // namespace ThermoFun

namespace Reaktoro {

// Forward declarations
class Element;
class AqueousSpecies;
class FluidSpecies;
using GaseousSpecies = FluidSpecies;
using LiquidSpecies = FluidSpecies;
class MineralSpecies;

/// Provides operations to retrieve physical and thermodynamic data of chemical species.
///
/// The Database class is used to retrieve information of chemical species. It is initialized
/// from a `xml` database file, which must satisfy some format conditions. Once it is initialized,
/// one can, for example, retrieve thermodynamic information of an aqueous species that will be
/// used to calculate its standard chemical potential.
///
//////*Usage**
///
/// In the example below, a Database instance is initialized and
/// queries are made to retrieve information of aqueous, gaseous
/// and liquid species.
/// Note that an exception is thrown if a species is not present in the database.
/// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// using namespace Reaktoro;
///
/// // Create a Database instance by parsing a local database file
/// Database database("geodb.xml")
///
/// // Retrieve information of species H2O(l), CO2(g) and CO2(liq)
/// AqueousSpecies aqueousSpecies = database.aqueousSpecies("H2O(l)");
/// GaseousSpecies gaseousSpecies = database.gaseousSpecies("CO2(g)");
/// LiquidSpecies liquidSpecies = database.liquidSpecies("CO2(liq)");

///
/// // Output the data of the species H2O(l), CO2(g), CO2(liq)
/// std::cout << aqueousSpecies << std::endl;
/// std::cout << gaseousSpecies << std::endl;
/// std::cout << liquidSpecies << std::endl;
/// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
///
/// @see AqueousSpecies, GaseousSpecies, LiquidSpecies, MineralSpecies
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
    /// exists with a given name, an exception will be thrown.
    /// @param filename The name of the database file
    explicit Database(std::string filename);

    Database(const ThermoFun::Database& fundatabase);

    /// Add an Element instance in the database.
    auto addElement(const Element& element) -> void;

    /// Add an AqueousSpecies instance in the database.
    auto addAqueousSpecies(const AqueousSpecies& species) -> void;

    /// Add a GaseousSpecies instance in the database.
    auto addGaseousSpecies(const GaseousSpecies& species) -> void;

    /// Add a LiquidSpecies instance in the database
    auto addLiquidSpecies(const LiquidSpecies& species) -> void;

    /// Add a MineralSpecies instance in the database.
    auto addMineralSpecies(const MineralSpecies& species) -> void;

    /// Return all elements in the database
    auto elements() const -> std::vector<Element>;

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

    /// Return all liquid species in the database
    auto liquidSpecies() -> std::vector<LiquidSpecies>;

    /// Return a liquid species in the database.
    /// **Note:** An exception is thrown if the database does not contain the species.
    /// @param name The name of the liquid species
    auto liquidSpecies(std::string name) const -> const LiquidSpecies&;

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

    /// Check if the database contains a given liquid species
    /// @param species The name of the liquid species
    auto containsLiquidSpecies(std::string species) const -> bool;

    /// Check if the database contains a given mineral species
    /// @param species The name of the mineral species
    auto containsMineralSpecies(std::string species) const -> bool;

    /// Return the aqueous species that contains at least one of the specified elements.
    auto aqueousSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<AqueousSpecies>;

    /// Return the gaseous species that contains at least one of the specified elements.
    auto gaseousSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<GaseousSpecies>;

    /// Return the liquid species that contains at least one of the specified elements.
    auto liquidSpeciesWithElements(const std::vector<std::string>& elements) const->std::vector<LiquidSpecies>;

    /// Return the mineral species that contains at least one of the specified elements.
    auto mineralSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<MineralSpecies>;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
