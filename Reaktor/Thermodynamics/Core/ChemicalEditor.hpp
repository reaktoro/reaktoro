// Reaktor is a C++ library for computational reaction modelling.
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

// C++ includes
#include <string>
#include <vector>
#include <memory>

// Reaktor includes

namespace Reaktor {

// Reaktor forward declarations
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
/// using namespace Reaktor;
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

    /// Add an aqueous phase in the chemical editor.
    /// Note that only one aqueous phase can exist in the chemical editor.
    /// So whenever this method is called, it has the effect of updating the
    /// current state of the aqueous phase in the editor.
    /// @param phase The AqueousPhase instance
    /// @return A reference to the aqueous phase
    auto addPhase(const AqueousPhase& phase) -> AqueousPhase&;

    /// Add a gaseous phase in the chemical editor.
    /// Note that only one gaseous phase can exist in the chemical editor.
    /// So whenever this method is called, it has the effect of updating the
    /// current state of the gaseous phase in the editor.
    /// @param phase The GaseousPhase instance
    /// @return A reference to the gaseous phase
    auto addPhase(const GaseousPhase& phase) -> GaseousPhase&;

    /// Add a mineral phase in the chemical editor.
    /// @param phase The MineralPhase instance
    /// @return A reference to the new mineral phase
    auto addPhase(const MineralPhase& phase) -> MineralPhase&;

    /// Add a mineral reaction in the chemical editor.
    /// @param The MineralReaction instance
    /// @return A reference to the new mineral reaction
    auto addReaction(const MineralReaction& reaction) -> MineralReaction&;

    /// Add an aqueous phase in the chemical editor.
    /// @param species The names of the species that compose the aqueous phase
    /// @return A reference to the aqueous phase
    auto addAqueousPhase(const std::vector<std::string>& species) -> AqueousPhase&;

    /// Add an aqueous phase in the chemical editor.
    /// @param species A string containing a list of species separated by space
    /// @return A reference to the aqueous phase
    auto addAqueousPhase(const std::string& species) -> AqueousPhase&;

    /// Add a gaseous phase in the chemical editor.
    /// @param species The names of the species that compose the gaseous phase
    /// @return A reference to the gaseous phase
    auto addGaseousPhase(const std::vector<std::string>& species) -> GaseousPhase&;

    /// Add a gaseous phase in the chemical editor.
    /// @param species A string containing a list of species separated by space
    /// @return A reference to the gaseous phase
    auto addGaseousPhase(const std::string& species) -> GaseousPhase&;

    /// Add a mineral phase in the chemical editor.
    /// @param species The names of the species that compose the mineral phase
    /// @return A reference to the new mineral phase
    auto addMineralPhase(const std::vector<std::string>& species) -> MineralPhase&;

    /// Add a mineral phase in the chemical editor.
    /// @param species A string containing a list of species separated by space
    /// @return A reference to the new mineral phase
    auto addMineralPhase(const std::string& species) -> MineralPhase&;

    /// Add a mineral reaction in the chemical editor.
    /// @param mineral The name of the mineral for which the reaction will be defined
    /// @return A reference to a MineralReaction instance
    auto addMineralReaction(const std::string& mineral) -> MineralReaction&;

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

} // namespace Reaktor
