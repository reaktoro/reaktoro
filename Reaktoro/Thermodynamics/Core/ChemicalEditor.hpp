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
#include <string>
#include <vector>
#include <memory>

// Reaktoro includes

namespace Reaktoro {

// Forward declarations
class Database;
class AqueousPhase;
class GaseousPhase;
class MineralPhase;
class ChemicalSystem;
class ReactionSystem;
class MineralReaction;

/// Provides convenient operations to initialize ChemicalSystem and ReactionSystem instances.
/// The ChemicalEditor class is used to conveniently create instances of classes ChemicalSystem and ReactionSystem.
///
/// **Usage**
///
/// The code below uses a ChemicalEditor instance to define a chemical system
/// composed of an aqueous phase and a mineral phase. In addition, it defines a
/// chemical reaction involving both aqueous and mineral species. Finally, the
/// editor instance is used to create a ChemicalSystem instance and a ReactionSystem
/// instance.
///
/// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// using namespace Reaktoro;
///
/// // A Database instance is necessary to create a ChemicalEditor instance
/// Database database("geodb.xml");
///
/// // Create the ChemicalEditor instance to setup the chemical system and reactions
/// ChemicalEditor editor(database);
///
/// // Define a chemical system with an aqueous and a mineral phase
/// editor.addAqueousPhase("H2O(l), H+, OH-, HCO3-, CO2(aq), Ca++");
/// editor.addMineralPhase("Calcite");
///
/// // Define a mineral reaction involving the mineral phase and the aqueous phase
/// editor.addMineralReaction("Calcite")
///     .setEquation("-1:Calcite, -1:H+, 1:HCO3-, 1:Ca++")
///     .setSpecificSurfaceArea(100.0, "cm2/g")
///     .addMechanism("logk = -5.81 mol/(m2*s), Ea = 23.5 kJ/mol")
///     .addMechanism("logk = -0.30 mol/(m2*s), Ea = 14.4 kJ/mol, a[H+] = 1.0");
///
/// // Create the ChemicalSystem and ReactionSystem instances
/// ChemicalSystem system = editor;
/// ReactionSystem reactions = editor;
/// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
///
/// @see Database, ChemicalSystem, ReactionSystem, AqueousPhase, GaseousPhase,
/// MineralPhase, AqueousSpecies, GaseousSpecies, MineralSpecies, MineralReaction
///
/// @ingroup Core
class ChemicalEditor
{
public:
    /// Construct a default ChemicalEditor instance.
    /// The built-in database `supcrt98.xml` is used for the initialization
    /// of a ChemicalEditor instance using the default constructor.
    ChemicalEditor();

    /// Construct a ChemicalEditor instance with the provided database.
    explicit ChemicalEditor(const Database& database);

    /// Construct a copy of the provided ChemicalEditor instance.
    ChemicalEditor(const ChemicalEditor& other);

    /// Destroy the ChemicalEditor instance.
    virtual ~ChemicalEditor();

    /// Assign this ChemicalEditor instance with another.
    auto operator=(const ChemicalEditor& other) -> ChemicalEditor&;

    /// Set the temperatures for constructing interpolation tables of thermodynamic properties.
    /// @param values The temperature values
    /// @param units The units of the temperature values
    auto setTemperatures(std::vector<double> values, std::string units) -> void;

    /// Set the pressures for constructing interpolation tables of thermodynamic properties.
    /// @param values The pressure values
    /// @param units The units of the pressure values
    auto setPressures(std::vector<double> values, std::string units) -> void;

    /// Initialize all possible phases that can exist with given elements.
    /// @param elements The element symbols of interest.
    auto initializePhasesWithElements(std::vector<std::string> elements) -> void;

    /// Add an aqueous phase in the chemical editor.
    /// Note that only one aqueous phase can exist in the chemical editor.
    /// So whenever this method is called, it has the effect of updating the
    /// current state of the aqueous phase in the editor.
    /// @param phase The AqueousPhase instance
    /// @return A reference to the created AqueousPhase object.
    auto addPhase(const AqueousPhase& phase) -> AqueousPhase&;

    /// Add a gaseous phase in the chemical editor.
    /// Note that only one gaseous phase can exist in the chemical editor.
    /// So whenever this method is called, it has the effect of updating the
    /// current state of the gaseous phase in the editor.
    /// @param phase The GaseousPhase instance
    /// @return A reference to the created GaseousPhase object.
    auto addPhase(const GaseousPhase& phase) -> GaseousPhase&;

    /// Add a mineral phase in the chemical editor.
    /// If a mineral phase with same name already exists, then
    /// the existing phase is replaced by the new one.
    /// @param phase The MineralPhase instance
    /// @return A reference to the created MineralPhase object.
    auto addPhase(const MineralPhase& phase) -> MineralPhase&;

    /// Add a mineral reaction in the chemical editor.
    /// @param reaction The MineralReaction instance
    /// @return A reference to the created MineralReaction object.
    auto addReaction(const MineralReaction& reaction) -> MineralReaction&;

    /// Add an aqueous phase in the chemical editor.
    /// This method constructs an AqueousPhase object that represents an aqueous phase in the system.
    /// The AqueousPhase object is created by specifying the names of the species one by one. These
    /// species names must conform to those used in the database that was specified during the
    /// initialization of the ChemicalEditor object, otherwise an exception will be thrown.
    /// The example below describes an usage of this method for an aqueous phase that could be
    /// formed by mixing H2O, CO2 and NaCl.
    /// ~~~
    /// #include <Reaktoro/Reaktoro.hpp> {delete}
    /// using namespace Reaktoro; {delete}
    /// ChemicalEditor editor;
    /// editor.addAqueousPhase({"H2O(l)", "H+", "OH-", "HCO3-", "CO2(aq)", "Na+", "Cl-"});
    /// ~~~
    /// An alternative way, in which to prior knowledge of the species names in the database is
    /// needed, consists of specifying a list of chemical element, compound, or substance names,
    /// and let the ChemicalEditor to figure out automatically which species from the loaded database
    /// should be added in the phase. This functionality is supported by method
    /// @ref addAqueousPhase(std::string compounds).
    /// @param species A vector containing the names of the species.
    /// @return A reference to the created AqueousPhase object.
    /// @see addGaseousPhase, addMineralPhase
    auto addAqueousPhase(std::vector<std::string> species) -> AqueousPhase&;

    /// @overload
    auto addAqueousPhase(std::initializer_list<std::string> species) -> AqueousPhase&;

    /// Add an aqueous phase in the chemical editor.
    /// This method constructs an AqueousPhase object that represents an aqueous phase in the system.
    /// Instead of listing the names of the species one by one, which might require prior knowledge
    /// of the species names in the database, this method permits the AqueousPhase object to be
    /// constructed by using a list of chemical element names or a list of compound or substance
    /// names that might not represent names of species in the database. The list of compounds
    /// will be broken into a list of element names, and the database will then be searched for all
    /// species that could be formed out of those elements. These species will then be used to
    /// construct the AqueousPhase object.
    /// The example below describes three equivalent alternatives to construct an AqueousPhase
    /// object that represents an aqueous phase that could be formed by mixing H2O, CO2 and NaCl.
    /// ~~~
    /// #include <Reaktoro/Reaktoro.hpp> {delete}
    /// using namespace Reaktoro; {delete}
    /// ChemicalEditor editor;
    /// editor.addAqueousPhase("H2O CO2 NaCl");
    /// editor.addAqueousPhase("HCNaCl");
    /// editor.addAqueousPhase("H C Na Cl");
    /// ~~~
    /// @note If only one name is given in `compounds`, and this name corresponds to a species
    /// in the database, then the phase will be created with only that species. The example below
    /// will produce an AqueousPhase with only species H2O(l), and not a phase with all possible
    /// species that could result from the combination of chemical elements H and O.
    /// ~~~
    /// #include <Reaktoro/Reaktoro.hpp> {delete}
    /// using namespace Reaktoro; {delete}
    /// ChemicalEditor editor;
    /// editor.addAqueousPhase("H2O(l)");
    /// ~~~
    /// This might not be relevant for an aqueous phase, which in general contains many species,
    /// but it is a convenient functionality for gaseous and mineral phases, for example, which
    /// might only contain one species of interest.
    /// @param compounds A string containing a list of element or compound names.
    /// @return A reference to the created AqueousPhase object.
    /// @see addGaseousPhase, addMineralPhase
    auto addAqueousPhase(std::string compounds) -> AqueousPhase&;

    /// Add a gaseous phase in the chemical editor.
    /// This method constructs a GaseousPhase object that represents a gaseous phase in the system.
    /// The GaseousPhase object is created by specifying the names of the species one by one. These
    /// species names must conform to those used in the database that was specified during the
    /// initialization of the ChemicalEditor object, otherwise an exception will be thrown.
    /// The example below describes an usage of this method for a gaseous phase that could be
    /// formed by mixing CH4 and O2.
    /// ~~~
    /// #include <Reaktoro/Reaktoro.hpp> {delete}
    /// using namespace Reaktoro; {delete}
    /// ChemicalEditor editor;
    /// editor.addGaseousPhase({"H2O(g)", "CO2(g)", "O2(g)", "CH4(g)"});
    /// ~~~
    /// An alternative way, in which to prior knowledge of the species names in the database is
    /// needed, consists of specifying a list of chemical element, compound, or substance names,
    /// and let the ChemicalEditor to figure out automatically which species from the loaded database
    /// should be added in the phase. This functionality is supported by method
    /// @ref addGaseousPhase(std::string compounds).
    /// @param species A vector containing the names of the species.
    /// @return A reference to the created GaseousPhase object.
    /// @see addAqueousPhase, addMineralPhase
    auto addGaseousPhase(std::vector<std::string> species) -> GaseousPhase&;

    /// @overload
    auto addGaseousPhase(std::initializer_list<std::string> species) -> GaseousPhase&;

    /// Add a gaseous phase in the chemical editor.
    /// This method constructs a GaseousPhase object that represents a gaseous phase in the system.
    /// Instead of listing the names of the species one by one, which might require prior knowledge
    /// of the species names in the database, this method permits the GaseousPhase object to be
    /// constructed by using a list of chemical element names or a list of compound or substance
    /// names that might not represent names of species in the database. The list of compounds
    /// will be broken into a list of element names, and the database will then be searched for all
    /// species that could be formed out of those elements. These species will then be used to
    /// construct the GaseousPhase object.
    /// The example below describes three equivalent alternatives to construct a GaseousPhase
    /// object that represents a gaseous phase that could be formed by mixing H2O and CO2.
    /// ~~~
    /// #include <Reaktoro/Reaktoro.hpp> {delete}
    /// using namespace Reaktoro; {delete}
    /// ChemicalEditor editor;
    /// editor.addGaseousPhase("H2O CO2");
    /// editor.addGaseousPhase("HOC");
    /// editor.addGaseousPhase("H O C");
    /// ~~~
    /// @note If only one name is given in `compounds`, and this name corresponds to a species
    /// in the database, then the phase will be created with only that species. The example below
    /// will produce a GaseousPhase with only species CO2(g), and not a phase with all possible
    /// species that could result from the combination of chemical elements C and O, such as CO(g),
    /// which could be achieved by specifying CO2 instead of CO2(g).
    /// ~~~
    /// #include <Reaktoro/Reaktoro.hpp> {delete}
    /// using namespace Reaktoro; {delete}
    /// ChemicalEditor editor;
    /// editor.addGaseousPhase("CO2(g)");
    /// ~~~
    /// @param compounds A string containing a list of element or compound names.
    /// @return A reference to the created GaseousPhase object.
    /// @see addAqueousPhase, addMineralPhase
    auto addGaseousPhase(std::string compounds) -> GaseousPhase&;

    /// Add a mineral phase in the chemical editor.
    /// This method constructs a MineralPhase object that represents a mineral phase in the system.
    /// The MineralPhase object is created by specifying the names of the species one by one. These
    /// species names must conform to those used in the database that was specified during the
    /// initialization of the ChemicalEditor object, otherwise an exception will be thrown.
    /// The example below describes an usage of this method for the creation of three pure mineral
    /// phases and one solid solution with two mineral species.
    /// ~~~
    /// #include <Reaktoro/Reaktoro.hpp> {delete}
    /// using namespace Reaktoro; {delete}
    /// ChemicalEditor editor;
    ///
    /// // Create a pure mineral phase with only calcite [CaCO3(s)]
    /// editor.addMineralPhase({"Calcite"});
    ///
    /// // Create a pure mineral phase with only magnesite [MgCO3(s)]
    /// editor.addMineralPhase({"Magnesite"});
    ///
    /// // Create a pure mineral phase with only dolomite [CaMg(CO3)2(s)]
    /// editor.addMineralPhase({"Dolomite"});
    ///
    /// // Create a solid solution with mineral species calcite and magnesite
    /// editor.addMineralPhase({"Calcite", "Magnesite"});
    /// ~~~
    /// An alternative way, in which to prior knowledge of the species names in the database is
    /// needed, consists of specifying a list of chemical element, compound, or substance names,
    /// and let the ChemicalEditor to figure out automatically which species from the loaded database
    /// should be added in the phase. This functionality is supported by method
    /// @ref addMineralPhase(std::string compounds).
    /// @param species A vector containing the names of the species.
    /// @return A reference to the created MineralPhase object.
    /// @see addAqueousPhase, addGaseousPhase
    auto addMineralPhase(std::vector<std::string> species) -> MineralPhase&;

    /// @overload
    auto addMineralPhase(std::initializer_list<std::string> species) -> MineralPhase&;

    /// Add a mineral phase in the chemical editor.
    /// This method constructs a MineralPhase object that represents a mineral phase in the system.
    /// Instead of listing the names of the species one by one, which might require prior knowledge
    /// of the species names in the database, this method permits the MineralPhase object to be
    /// constructed by using a list of chemical element names or a list of compound or substance
    /// names that might not represent names of species in the database. The list of compounds
    /// will be broken into a list of element names, and the database will then be searched for all
    /// species that could be formed out of those elements. These species will then be used to
    /// construct the MineralPhase object.
    /// The example below describes several possibilities to construct a MineralPhase object.
    /// ~~~
    /// #include <Reaktoro/Reaktoro.hpp> {delete}
    /// using namespace Reaktoro; {delete}
    /// ChemicalEditor editor;
    ///
    /// // Create a solid solution with all minerals that could be formed by combining compounds CaCO3 and MgCO3.
    /// editor.addMineralPhase("CaCO3 MgCO3");
    ///
    /// // Create a pure mineral phase with only calcite [CaCO3(s)]
    /// editor.addMineralPhase("Calcite"); // assuming the name Calcite is in the database
    ///
    /// // Create a solid solution with all minerals that could be formed from elements Ca, C, and O.
    /// editor.addMineralPhase("CaCO3");  // assuming the name CaCO3 is not in the database
    /// editor.addMineralPhase("Ca C O"); // equivalent to the previous call
    /// ~~~
    /// @note In most cases, the solid solutions of interest have predefined mineral composition, so that
    /// one might prefer instead to list the mineral end-members one by one, instead of letting
    /// ChemicalEditor to populate the solid solution with many minerals that could be formed from a
    /// given list of chemical elements or compound names.
    /// @param compounds A string containing a list of element or compound names.
    /// @return A reference to the created MineralPhase object.
    /// @see addAqueousPhase, addGaseousPhase
    auto addMineralPhase(std::string compounds) -> MineralPhase&;

    /// Add a mineral reaction in the chemical editor.
    /// @param reaction The mineral reaction.
    /// @return A reference to the created MineralReaction object.
    auto addMineralReaction(const MineralReaction& reaction) -> MineralReaction&;

    /// Add a mineral reaction in the chemical editor.
    /// @param mineral The name of the mineral for which the reaction will be defined.
    /// @return A reference to the created MineralReaction object.
    auto addMineralReaction(std::string mineral) -> MineralReaction&;

    /// Return the aqueous phase in the chemical editor.
    auto aqueousPhase() const -> const AqueousPhase&;

    /// Return the aqueous phase in the chemical editor.
    auto aqueousPhase() -> AqueousPhase&;

    /// Return the gaseous phase in the chemical editor.
    auto gaseousPhase() const -> const GaseousPhase&;

    /// Return the gaseous phase in the chemical editor.
    auto gaseousPhase() -> GaseousPhase&;

    /// Return the mineral phases in the chemical editor.
    auto mineralPhases() const -> const std::vector<MineralPhase>&;

    /// Return the mineral phases in the chemical editor.
    auto mineralPhases() -> std::vector<MineralPhase>&;

    /// Create a ChemicalSystem instance with the current state of the chemical editor
    auto createChemicalSystem() const -> ChemicalSystem;

    /// Create a ReactionSystem instance with the current state of the chemical editor
    auto createReactionSystem() const -> ReactionSystem;

    /// Convert this ChemicalEditor instance to a ChemicalSystem instance
    operator ChemicalSystem() const;

    /// Convert this ChemicalEditor instance to a ReactionSystem instance
    operator ReactionSystem() const;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
