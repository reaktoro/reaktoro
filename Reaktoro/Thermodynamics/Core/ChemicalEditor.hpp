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
class AqueousPhase;
class ChemicalSystem;
class Database;
class GaseousPhase;
class LiquidPhase;
class MineralPhase;
class MineralReaction;
class ReactionSystem;
class StringList;

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
/// LiquidPhase, MineralPhase, AqueousSpecies, GaseousSpecies, LiquidSpecies,
/// MineralSpecies, MineralReaction
///
/// @ingroup Core
class ChemicalEditor
{
public:
    /// Construct a default ChemicalEditor instance.
    /// The built-in database `supcrt98.xml` is used for the initialization
    /// of a ChemicalEditor instance using the default constructor.
    ChemicalEditor();

    /// Construct a ChemicalEditor instance with a provided database.
    explicit ChemicalEditor(const Database& database);

    /// Construct a ChemicalEditor instance with a provided ThermoFun database.
    explicit ChemicalEditor(const ThermoFun::Database& database);

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
    auto initializePhasesWithElements(const StringList& elements) -> void;

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
    /// current state on gaseous phase in the editor.
    /// @param phase The GaseousPhase instance
    /// @return A reference to the created GaseousPhase object.
    auto addPhase(const GaseousPhase& phase) -> GaseousPhase&;

    /// Add a liquid phase in the chemical editor.
    /// Note that only one liquid phase can exist in the chemical editor.
    /// So whenever this method is called, it has the effect of updating the
    /// current state on liquid phase in the editor.
    /// @param phase The LiquidPhase instance
    /// @return A reference to the created LiquidPhase object.
    auto addPhase(const LiquidPhase& phase) -> LiquidPhase&;

    /// Add a mineral phase in the chemical editor.
    /// If a mineral phase with the same name already exists, then
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
    /// initialization of the ChemicalEditor object, otherwise, an exception will be thrown.
    /// The example below describes the usage of this method for an aqueous phase that could be
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
    /// should be added in the phase. This functionality is supported by methods
    ///
    /// @ref addAqueousPhaseWithElements(const std::vector<std::string>& elements),
    ///
    /// @ref addAqueousPhaseWithElementsOf(const std::vector<std::string>& compounds).
    ///
    /// @param species A StringList containing the names of the species.
    /// @return A reference to the created AqueousPhase object.
    /// @see addGaseousPhase, addLiquidPhase, addMineralPhase
    ///
    /// @note The old use of this function to add elements and/or compounds was removed. To use these
    /// functionalities, use addAqueousPhaseWitElements to add elements and addAqueousPhaseWitElementsOf
    /// to add compounds.
    auto addAqueousPhase(const StringList& species) -> AqueousPhase&;

    /// Add an aqueous phase in the chemical editor.
    /// This method constructs an AqueousPhase object that represents an aqueous phase in the system.
    /// Instead of listing the names of the species one by one, which might require prior knowledge
    /// of the species names in the database, this method permits the AqueousPhase object to be
    /// constructed by using a list of chemical element names, and the database will then be searched for all
    /// species that could be formed out of those elements. These species will then be used to
    /// construct the AqueousPhase object.
    /// The example below describes three equivalent alternatives to construct an AqueousPhase
    /// object that represents an aqueous phase that could be formed by mixing H2O, CO2 and NaCl.
    /// ~~~
    /// ChemicalEditor editor;
    /// editor.addAqueousPhaseWithElements({"H", "C", "Na", "Cl"});
    /// editor.addAqueousPhaseWithElements({"H2O", "CO2", "NaCl");
    /// editor.addAqueousPhaseWithElements("HCNaCl");
    /// ~~~
    /// This might not be relevant for an aqueous phase, which in general contains many species,
    /// but it is a convenient functionality for gaseous and mineral phases, for example, which
    /// might only contain one species of interest.
    /// @param elements A StringList containing a list of chemical element names.
    /// @return A reference to the created AqueousPhase object.
    /// @see addGaseousPhaseWithElements, addLiquidPhaseWithElements, addMineralPhaseWithElements
    auto addAqueousPhaseWithElements(const StringList& elements) -> AqueousPhase&;

    /// Add an aqueous phase in the chemical editor.
    /// This method constructs an AqueousPhase object that represents an aqueous phase in the system.
    /// Instead of listing the names of the species one by one, which might require prior knowledge
    /// of the species names in the database, this method permits the AqueousPhase object to be
    /// constructed by using a list of compound or substance names that might not represent names of
    /// species in the database. The list of compounds will be broken into a list of element names,
    /// and the database will then be searched for all species that could be formed out of those elements.
    /// These species will then be used to construct the AqueousPhase object.
    /// The example below describes three equivalent alternatives to construct an AqueousPhase
    /// object that represents an aqueous phase that could be formed by mixing H2O, CO2 and NaCl.
    /// ~~~
    /// ChemicalEditor editor;
    /// editor.addAqueousPhaseWithElementsOf("H2O NaCl CO2");
    /// editor.addAqueousPhaseWithElementsOf({"H2O", "CO2", "NaCl");
    /// editor.addAqueousPhaseWithElementsOf("HCNaCl");
    /// ~~~
    /// This might not be relevant for an aqueous phase, which in general contains many species,
    /// but it is a convenient functionality for gaseous and mineral phases, for example, which
    /// might only contain one species of interest.
    /// @param compounds A StringList containing a list of compound names.
    /// @return A reference to the created AqueousPhase object.
    /// @see addGaseousPhaseWithElementsOf, addMineralPhaseWithElementsOf
    auto addAqueousPhaseWithElementsOf(const StringList& compounds)->AqueousPhase&;

    /// Add a gaseous phase in the chemical editor.
    /// This method constructs a GaseousPhase object that represents a gaseous phase in the system.
    /// The GaseousPhase object is created by specifying the names of the species one by one. These
    /// species names must conform to those used in the database that was specified during the
    /// initialization of the ChemicalEditor object, otherwise, an exception will be thrown.
    /// The example below describes the usage of this method for a gaseous phase that could be
    /// formed by mixing CH4 and O2.
    /// ~~~
    /// ChemicalEditor editor;
    /// editor.addGaseousPhase({"H2O(g)", "CO2(g)", "O2(g)", "CH4(g)"});
    /// ~~~
    /// An alternative way, in which to prior knowledge of the species names in the database is
    /// needed, consists of specifying a list of chemical element, compound, or substance names,
    /// and let the ChemicalEditor to figure out automatically which species from the loaded database
    /// should be added in the phase. This functionality is supported by method
    /// @ref addGaseousPhaseWithElements(std::string elements).
    /// @ref addGaseousPhaseWithElementsOf(std::string compounds).
    /// @param species A StringList containing the names of the species.
    /// @return A reference to the created GaseousPhase object.
    /// @see addAqueousPhase, addLiquidPhase, addMineralPhase
    ///
    /// @note The old use of this function to add elements and/or compounds was removed. To use these
    /// functionalities, use addGaseousPhaseWitElements to add elements and addGaseousPhaseWitElementsOf
    /// to add compounds.
    auto addGaseousPhase(const StringList& species) -> GaseousPhase&;

    /// Add a gaseous phase in the chemical editor.
    /// This method constructs a GaseousPhase object that represents a gaseous phase in the system.
    /// Instead of listing the names of the species one by one, which might require prior knowledge
    /// of the species names in the database, this method permits the GaseousPhase object to be
    /// constructed by using a list of chemical element names and the database will then be searched for all
    /// species that could be formed out of those elements. These species will then be used to
    /// construct the GaseousPhase object.
    /// The example below describes three equivalent alternatives to construct a GaseousPhase
    /// object that represents a gaseous phase that could be formed by mixing H2O and CO2.
    /// ~~~
    /// ChemicalEditor editor;
    /// editor.addGaseousPhaseWithElements({"H", "O", "C"});
    /// editor.addGaseousPhaseWithElements({"H2O", "CO2"});
    /// editor.addGaseousPhaseWithElements({"HOC"});
    /// ~~~
    /// @param elements A StringList containing a list of chemical element names.
    /// @return A reference to the created GaseousPhase object.
    /// @see addAqueousPhaseWithElements, addLiquidPhaseWithElements, addMineralPhaseWithElements
    auto addGaseousPhaseWithElements(const StringList& elements) -> GaseousPhase&;

    /// Add a gaseous phase in the chemical editor.
    /// This method constructs a GaseousPhase object that represents a gaseous phase in the system.
    /// Instead of listing the names of the species one by one, which might require prior knowledge
    /// of the species names in the database, this method permits the GaseousPhase object to be
    /// constructed by using a list of compound or substance names that might not represent names of
    /// species in the database. The list of compounds will be broken into a list of element names,
    /// and the database will then be searched for all species that could be formed out of those elements.
    // These species will then be used to construct the GaseousPhase object.
    /// The example below describes three equivalent alternatives to construct a GaseousPhase
    /// object that represents a gaseous phase that could be formed by mixing H2O and CO2.
    /// ~~~
    /// ChemicalEditor editor;
    /// editor.addGaseousPhaseWithElementsOf({"H2O", "CO2"});
    /// editor.addGaseousPhaseWithElementsOf({"H2O CO2"});
    /// editor.addGaseousPhaseWithElementsOf({"HOC"});
    /// ~~~
    /// @param compounds A StringList containing a list of compound names.
    /// @return A reference to the created GaseousPhase object.
    /// @see addAqueousPhaseWithElements, addLiquidPhaseWithElementsOf, addMineralPhaseWithElements
    auto addGaseousPhaseWithElementsOf(const StringList& compounds) -> GaseousPhase&;

    /// Add a liquid phase in the chemical editor.
    /// This method constructs a LiquidPhase object that represents a liquid phase in the system.
    /// The LiquidPhase object is created by specifying the names of the species one by one. These
    /// species names must conform to those used in the database that was specified during the
    /// initialization of the ChemicalEditor object, otherwise, an exception will be thrown.
    /// The example below describes the usage of this method for a liquid phase that could be
    /// formed by mixing CH4 and O2.
    /// ~~~
    /// ChemicalEditor editor;
    /// editor.addLiquidPhase({"H2O(liq)", "CO2(liq)", "O2(liq)", "CH4(liq)"});
    /// ~~~
    /// An alternative way, in which to prior knowledge of the species names in the database is
    /// needed, consists of specifying a list of chemical element, compound, or substance names,
    /// and let the ChemicalEditor to figure out automatically which species from the loaded database
    /// should be added in the phase. This functionality is supported by method
    /// @ref addLiquidPhaseWithElements(std::string elements).
    /// @ref addLiquidPhaseWithElementsOf(std::string compounds).
    /// @param species A StringList containing the names of the species.
    /// @return A reference to the created LiquidPhase object.
    /// @see addAqueousPhase, addGaseousPhase, addMineralPhase
    ///
    /// @note The old use of this function to add elements and/or compounds was removed. To use these
    /// functionalities, use addGaseousPhaseWitElements to add elements and addGaseousPhaseWitElementsOf
    /// to add compounds.
    auto addLiquidPhase(const StringList& species) -> LiquidPhase&;

    /// Add a liquid phase in the chemical editor.
    /// This method constructs a LiquidPhase object that represents a liquid phase in the system.
    /// Instead of listing the names of the species one by one, which might require prior knowledge
    /// of the species names in the database, this method permits the LiquidPhase object to be
    /// constructed by using a list of chemical element names and the database will then be searched for all
    /// species that could be formed out of those elements. These species will then be used to
    /// construct the LiquidPhase object.
    /// The example below describes three equivalent alternatives to construct a LiquidPhase
    /// object that represents a liquid phase that could be formed by mixing H2S and CO2.
    /// ~~~
    /// ChemicalEditor editor;
    /// editor.addLiquidPhaseWithElements({"H", "O", "C", "S"});
    /// editor.addLiquidPhaseWithElements({"H2S", "CO2"});
    /// editor.addLiquidPhaseWithElements({"HOCS"});
    /// ~~~
    /// @param elements A StringList containing a list of chemical element names.
    /// @return A reference to the created LiquidPhase object.
    /// @see addAqueousPhaseWithElements, addGaseousPhaseWithElements, addMineralPhaseWithElements
    auto addLiquidPhaseWithElements(const StringList& elements) -> LiquidPhase&;

    /// Add a liquid phase in the chemical editor.
    /// This method constructs a LiquidPhase object that represents a liquid phase in the system.
    /// Instead of listing the names of the species one by one, which might require prior knowledge
    /// of the species names in the database, this method permits the LiquidPhase object to be
    /// constructed by using a list of compound or substance names that might not represent names of
    /// species in the database. The list of compounds will be broken into a list of element names,
    /// and the database will then be searched for all species that could be formed out of those elements.
    // These species will then be used to construct the LiquidPhase object.
    /// The example below describes three equivalent alternatives to construct a LiquidPhase
    /// object that represents a liquid phase that could be formed by mixing H2S and CO2.
    /// ~~~
    /// ChemicalEditor editor;
    /// editor.addLiquidPhaseWithElementsOf({"H", "O", "C", "S"});
    /// editor.addLiquidPhaseWithElementsOf({"H2S", "CO2"});
    /// editor.addLiquidPhaseWithElementsOf({"HOCS"});
    /// ~~~
    /// @param compounds A StringList containing a list of compound names.
    /// @return A reference to the created LiquidPhase object.
    /// @see addAqueousPhaseWithElements, addGaseousPhaseWithElementsOf, addMineralPhaseWithElements
    auto addLiquidPhaseWithElementsOf(const StringList& compounds) -> LiquidPhase&;

    /// Add a mineral phase in the chemical editor.
    /// This method constructs a MineralPhase object that represents a mineral phase in the system.
    /// The MineralPhase object is created by specifying the names of the species one by one. These
    /// species names must conform to those used in the database that was specified during the
    /// initialization of the ChemicalEditor object, otherwise, an exception will be thrown.
    /// The example below describes the usage of this method for the creation of two pure mineral
    /// phases and one solid solution with two mineral species.
    /// ~~~
    /// ChemicalEditor editor;
    ///
    /// // Create a pure mineral phase with only calcite [CaCO3(s)]
    /// editor.addMineralPhase({"Calcite"});
    ///
    /// // Create a pure mineral phase with only magnesite [MgCO3(s)]
    /// editor.addMineralPhase("Magnesite");
    ///
    /// // Create a solid solution with mineral species calcite and magnesite
    /// editor.addMineralPhase({"Calcite", "Magnesite"});
    /// ~~~
    /// An alternative way, in which a prior knowledge of the species names in the database is
    /// needed, consists of specifying a list of chemical elements or compounds,
    /// and let the ChemicalEditor to figure out automatically which species from the loaded database
    /// should be added in the phase. This functionality is supported by methods
    /// @ref addMineralPhaseWithElements(std::string elements) and
	/// @ref addMineralPhaseWithElementsOf(std::string compounds).
    /// @param species A StringList containing the names of the species.
    /// @return A reference to the created MineralPhase object.
    /// @see addAqueousPhase, addLiquidPhase, addGaseousPhase
    ///
    /// @note The old use of this function to add elements and/or compounds was removed. To use these
    /// functionalities, use addMineralPhaseWitElements to add eslements and addMineralPhaseWitElementsOf
    /// to add compounds.
    auto addMineralPhase(const StringList& species) -> MineralPhase&;

    /// Add a mineral phase in the chemical editor.
    /// This method constructs a MineralPhase object that represents a mineral phase in the system.
    /// Instead of listing the names of the species one by one, which might require prior knowledge
    /// of the species names in the database, this method permits the MineralPhase object to be
    /// constructed by using a list of chemical element names, and the database will then be searched for all
    /// species that could be formed out of those elements. These species will then be used to
    /// construct the MineralPhase object.
    /// The example below describes several possibilities to construct a MineralPhase object.
    /// ~~~
    /// ChemicalEditor editor;
    ///
    /// // This will return an error, as it only accepts element names
    /// editor.addMineralPhaseWithElements({"CaCO3", "MgCO3"});
    ///
    /// editor.addMineralPhaseWithElements({"Ca", "C", "O"});
    ///
    /// // This will only recognize the element "O", and CaC will be ignored
    /// editor.addMineralPhaseWithElements({"CaC", "O"});
    /// ~~~
    /// @note In most cases, the solid solutions of interest have predefined mineral composition, so that
    /// one might prefer instead to list the mineral end-members one by one, instead of letting
    /// ChemicalEditor to populate the solid solution with many minerals that could be formed from a
    /// given list of chemical elements.
    /// @param elements A StringList containing a list of chemical element names.
    /// @return A reference to the created MineralPhase object.
    /// @see addAqueousPhaseWithElements, addLiquidPhaseWithElements, addGaseousPhaseWithElements
    auto addMineralPhaseWithElements(const StringList& elements) -> MineralPhase&;

    /// Add a mineral phase in the chemical editor.
    /// This method constructs a MineralPhase object that represents a mineral phase in the system.
    /// Instead of listing the names of the species one by one, which might require prior knowledge
    /// of the species names in the database, this method permits the MineralPhase object to be
    /// constructed by using a list of compound or substance names that might not represent names of
    /// species in the database. The list of compounds will be broken into a list of element names,
    /// and the database will then be searched for all species that could be formed out of those elements.
    /// These species will then be used to construct the MineralPhase object.
    /// The example below describes several possibilities to construct a MineralPhase object.
    /// ~~~
    /// ChemicalEditor editor;
    ///
    /// // Create a solid solution with all minerals that could be formed by combining compounds CaCO3 and MgCO3.
    /// editor.addMineralPhaseWithElementsOf({"CaCO3", "MgCO3"});
    ///
    /// // Create a solid solution with all minerals that could be formed from elements Ca, C, and O.
    /// editor.addMineralPhaseWithElementsOf("CaCO3");  // assuming the name CaCO3 is not in the database
    /// ~~~
    /// @note In most cases, the solid solutions of interest have predefined mineral composition, so that
    /// one might prefer instead to list the mineral end-members one by one, instead of letting
    /// ChemicalEditor to populate the solid solution with many minerals that could be formed from a
    /// given list of compound names.
    /// @param compounds A StringList containing a list of compound names.
    /// @return A reference to the created MineralPhase object.
    /// @see addAqueousPhaseWithElementsOf, addGaseousPhaseWithElementsOf
    auto addMineralPhaseWithElementsOf(const StringList& compounds) -> MineralPhase&;

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
    auto gaseousPhase() ->GaseousPhase&;

    /// Return the liquid phase in the chemical editor.
    auto liquidPhase() const -> const LiquidPhase&;

    /// Return the liquid phase in the chemical editor.
    auto liquidPhase() ->LiquidPhase&;

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
