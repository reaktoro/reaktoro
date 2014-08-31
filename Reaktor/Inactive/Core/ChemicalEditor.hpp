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
#include <string>
#include <vector>
#include <memory>

// Reaktor includes
#include <Reaktor/Core/Types.hpp>

namespace Reaktor {

// Forward declarations
class Database;
class AqueousPhase;
class GaseousPhase;
class MineralPhase;
class AqueousSpecies;
class GaseousSpecies;
class MineralSpecies;
class ChemicalSystem;
class ReactionSystem;
class MineralReaction;

/**
 * Provides convenient operations to initialize ChemicalSystem and ReactionSystem instances
 *
 * The ChemicalEditor class is used to conveniently create instances of classes
 * ChemicalSystem and ReactionSystem.
 *
 * **Usage**
 *
 * The code below uses a ChemicalEditor instance to define a chemical system
 * composed of an aqueous phase and a mineral phase. In addition, it defines a
 * chemical reaction involving both aqueous and mineral species. Finally, the
 * editor instance is used to create a ChemicalSystem instance and a ReactionSystem
 * instance.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * using namespace Reaktor;
 *
 * // A Database instance is necessary to create a ChemicalEditor instance
 * Database database("geodb.xml");
 *
 * // Create the ChemicalEditor instance to setup the chemical system and reactions
 * ChemicalEditor editor(database);
 *
 * // Define a chemical system with an aqueous and a mineral phase
 * editor.addAqueousPhase("H2O(l), H+, OH-, HCO3-, CO2(aq), Ca++");
 * editor.addMineralPhase("Calcite");
 *
 * // Define a mineral reaction involving the mineral phase and the aqueous phase
 * editor.addMineralReaction("Calcite")
 *     .setEquation("-1:Calcite, -1:H+, 1:HCO3-, 1:Ca++")
 *     .setSpecificSurfaceArea(100.0, "cm2/g")
 *     .addMechanism("logk = -5.81 mol/(m2*s), Ea = 23.5 kJ/mol")
 *     .addMechanism("logk = -0.30 mol/(m2*s), Ea = 14.4 kJ/mol, a[H+] = 1.0");
 *
 * // Create the ChemicalSystem and ReactionSystem instances
 * ChemicalSystem system = editor;
 * ReactionSystem reactions = editor;
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * @see Database, ChemicalSystem, ReactionSystem, AqueousPhase, GaseousPhase,
 * MineralPhase, AqueousSpecies, GaseousSpecies, MineralSpecies, MineralReaction
 *
 * @ingroup Core
 */
class ChemicalEditor
{
public:
    /**
     * Constructs a chemical editor instance with the provided database
     */
    explicit ChemicalEditor(const Database& database);

    /**
     * Constructs a copy of the provided chemical editor instance
     */
    ChemicalEditor(const ChemicalEditor& other);

    /**
     * Destructs the chemical editor instance
     */
    virtual ~ChemicalEditor();

    /**
     * Assigns this chemical editor instance with another
     */
    auto operator=(const ChemicalEditor& other) -> ChemicalEditor&;

    /**
     * Adds an aqueous phase in the chemical editor
     *
     * Note that only one aqueous phase can exist in the chemical editor.
     * So whenever this method is called, it has the effect of updating the
     * current state of the aqueous phase in the editor.
     *
     * @param phase The @ref AqueousPhase instance
     * @return A reference to the aqueous phase
     */
    auto addPhase(const AqueousPhase& phase) -> AqueousPhase&;

    /**
     * Adds a gaseous phase in the chemical editor
     *
     * Note that only one gaseous phase can exist in the chemical editor.
     * So whenever this method is called, it has the effect of updating the
     * current state of the gaseous phase in the editor.
     *
     * @param phase The @ref GaseousPhase instance
     * @return A reference to the gaseous phase
     */
    auto addPhase(const GaseousPhase& phase) -> GaseousPhase&;

    /**
     * Adds a mineral phase in the chemical editor
     * @param phase The @ref MineralPhase instance
     * @return A reference to the new mineral phase
     */
    auto addPhase(const MineralPhase& phase) -> MineralPhase&;

    /**
     * Adds a mineral reaction in the chemical editor
     * @param The @ref MineralReaction instance
     * @return A reference to the new mineral reaction
     */
    auto addReaction(const MineralReaction& reaction) -> MineralReaction&;

    /**
     * Adds an aqueous phase in the chemical editor containing all aqueous species in the loaded database
     * @return A reference to the aqueous phase
     */
    auto addAqueousPhase() -> AqueousPhase&;

    /**
     * Adds an aqueous phase in the chemical editor
     * @param species The aqueous species that compose the aqueous phase
     * @return A reference to the aqueous phase
     */
    auto addAqueousPhase(const std::vector<AqueousSpecies>& species) -> AqueousPhase&;

	/**
     * Adds an aqueous phase in the chemical editor
     * @param species The names of the species that compose the aqueous phase
     * @return A reference to the aqueous phase
     */
    auto addAqueousPhase(const std::vector<std::string>& species) -> AqueousPhase&;

    /**
     * Adds an aqueous phase in the chemical editor
     * @param species A string containing a list of species separated by space
     * @return A reference to the aqueous phase
     */
    auto addAqueousPhase(const std::string& species) -> AqueousPhase&;

    /**
     * Adds a gaseous phase in the chemical editor containing all gaseous species in the loaded database
     * @return A reference to the gaseous phase
     */
    auto addGaseousPhase() -> GaseousPhase&;

    /**
     * Adds a gaseous phase in the chemical editor
     * @param species The gaseous species that compose the gaseous phase
     * @return A reference to the gaseous phase
     */
    auto addGaseousPhase(const std::vector<GaseousSpecies>& species) -> GaseousPhase&;

	/**
     * Adds a gaseous phase in the chemical editor
     * @param species The names of the species that compose the gaseous phase
     * @return A reference to the gaseous phase
     */
	auto addGaseousPhase(const std::vector<std::string>& species) -> GaseousPhase&;

	/**
     * Adds a gaseous phase in the chemical editor
     * @param species A string containing a list of species separated by space
     * @return A reference to the gaseous phase
     */
    auto addGaseousPhase(const std::string& species) -> GaseousPhase&;

	/**
     * Adds a mineral phase in the chemical editor
     * @param species The mineral species that compose the mineral phase
     * @return A reference to the new mineral phase
     */
	auto addMineralPhase(const std::vector<MineralSpecies>& species) -> MineralPhase&;

	/**
     * Adds a mineral phase in the chemical editor
     * @param species The names of the species that compose the mineral phase
     * @return A reference to the new mineral phase
     */
	auto addMineralPhase(const std::vector<std::string>& species) -> MineralPhase&;

	/**
     * Adds a mineral phase in the chemical editor
     * @param species A string containing a list of species separated by space
     * @return A reference to the new mineral phase
     */
    auto addMineralPhase(const std::string& species) -> MineralPhase&;

    /**
     * Adds a mineral reaction in the chemical editor
     * @param mineral The name of the mineral for which the reaction will be defined
     * @return A reference to a @ref MineralReaction instance
     */
    auto addMineralReaction(const std::string& mineral) -> MineralReaction&;

    /**
     * Sets the standard chemical potential function of the species
     *
     * Setting the standard chemical potential function of a species
     * using this method overwrites any other function already set.
     *
     * @param species The name of the chemical species
     * @param potential The standard chemical potential function
     */
    auto setChemicalPotentialFn(const std::string& species, const ChemicalPotentialFn& potential) -> void;

    /**
     * Gets the aqueous phase in the chemical editor
     */
    auto aqueousPhase() const -> const AqueousPhase&;

    /**
     * Gets the aqueous phase in the chemical editor
     */
    auto aqueousPhase() -> AqueousPhase&;

    /**
     * Gets the gaseous phase in the chemical editor
     */
    auto gaseousPhase() const -> const GaseousPhase&;

    /**
     * Gets the gaseous phase in the chemical editor
     */
    auto gaseousPhase() -> GaseousPhase&;

    /**
     * Gets the mineral phases in the chemical editor
     */
    auto mineralPhases() const -> const std::vector<MineralPhase>&;

    /**
     * Gets the mineral phases in the chemical editor
     */
    auto mineralPhases() -> std::vector<MineralPhase>&;

    /**
     * Creates a @ref ChemicalSystem instance with the current state of the chemical editor
     */
	auto createChemicalSystem() const -> ChemicalSystem;

	/**
	 * Creates a @ref ReactionSystem instance with the current state of the chemical editor
	 */
	auto createReactionSystem() const -> ReactionSystem;

	/**
     * Converts this chemical editor to a @ref ChemicalSystem instance
     */
	operator ChemicalSystem() const;

	/**
	 * Converts this chemical editor to a @ref ReactionSystem instance
	 */
	operator ReactionSystem() const;

private:
	/// The implementation class
	class Impl;

	/// The pointer to the implementation instance
	std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktor
