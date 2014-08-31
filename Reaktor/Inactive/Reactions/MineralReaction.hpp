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
#include <Reaktor/Common/Units.hpp>
#include <Reaktor/Core/ReactionEquation.hpp>
#include <Reaktor/Core/Types.hpp>
#include <Reaktor/Reactions/MineralCatalyst.hpp>
#include <Reaktor/Reactions/MineralMechanism.hpp>

namespace Reaktor {

// Forward declarations
class ChemicalSystem;
class Reaction;

class MineralReaction
{
public:
    /**
     * Constructs a default @ref MineralReaction instance
     */
    MineralReaction();

    /**
     * Constructs a @ref MineralReaction instance with the mineral name
     */
    MineralReaction(const std::string& mineral);

    /**
     * Constructs a copy of a MineralReaction instance
     */
    MineralReaction(const MineralReaction& other);

    /**
     * Destroys this instance
     */
    virtual ~MineralReaction();

    /**
     * Assigns a MineralReaction instance to this instance
     */
    auto operator=(MineralReaction other) -> MineralReaction&;

    /**
     * Sets the name of the mineral species in the reaction
     */
    auto setMineral(const std::string& mineral) -> MineralReaction&;

    /**
     * Sets the equation of the mineral reaction
     *
     * @see ReactionEquation
     */
    auto setEquation(const ReactionEquation& equation) -> MineralReaction&;

    /**
     * Sets the equation of the mineral reaction
     *
     * The format of @c equation can be found in @ref ReactionEquation.
     *
     * @param equation The @c string defining the reaction equation
     *
     * @see ReactionEquation
     */
    auto setEquation(const std::string& equation) -> MineralReaction&;

    /**
     * Sets the equilibrium constant of the mineral reaction
     *
     * If no equilibrium contant is provided, it will be calculated from the chemical potentials
     * of the participating species in the reaction.
     *
     * @see EquilibriumConstant
     */
    auto setEquilibriumConstant(const EquilibriumConstantFn& equilibrium_constant) -> MineralReaction&;

    /**
     * Sets the specific surface area of the mineral
     *
     * The specific surface area of the mineral can be set using units
     * that are convertible to either m<sup>2</sup>/g or m<sup>2</sup>/m<sup>3</sup>.
     * Whichever unit is used, the specified specific surface area will be automatically
     * converted to a molar surface area with units of m<sup>2</sup>/mol.
     *
     * @param value The value of the specific surface area
     * @param unit The units of the specific surface area (must be convertible to either m2/g or m2/m3)
     */
    auto setSpecificSurfaceArea(double value, const std::string& unit) -> MineralReaction&;

    /**
     * Adds a mineral mechanism to the kinetic rate model of the mineral reaction
     * @see MineralMechanism
     */
    auto addMechanism(const std::string& mechanism) -> MineralReaction&;

    /**
     * Adds a mineral mechanism to the kinetic rate model of the mineral reaction
     *
     * @param mechanism The mechanism to be considered in the kinetic rate model
     *
     * @see MineralMechanism
     */
    auto addMechanism(const MineralMechanism& mechanism) -> MineralReaction&;

    /**
     * Sets the mineral mechanisms of the kinetic rate model of the mineral reaction
     *
     * @param mechanisms The mechanisms of the kinetic rate model
     *
     * @see MineralMechanism
     */
    auto setMechanisms(const std::vector<MineralMechanism>& mechanisms) -> MineralReaction&;

    /**
     * Gets the name of the mineral species in the reaction
     */
    auto mineral() const -> const std::string&;

    /**
     * Gets the equation of the mineral reaction as a @c string
     *
     * @see ReactionEquation
     */
    auto equation() const -> const ReactionEquation&;

    /**
     * Gets the equilibrium constant of the mineral reaction
     *
     * @see EquilibriumConstant
     */
    auto equilibriumConstant() const -> const EquilibriumConstantFn&;

    /**
     * Gets the specific surface area of the mineral
     */
    auto specificSurfaceArea() const -> units::SpecificSurfaceArea;

    /**
     * Gets the volumetric surface area of the mineral
     */
    auto volumetricSurfaceArea() const -> units::VolumetricSurfaceArea;

    /**
     * Gets the mineral mechanisms of the kinetic rate model of the mineral reaction
     *
     * @see MineralMechanism
     */
    auto mechanisms() const -> const std::vector<MineralMechanism>&;

private:
    class Impl;

    std::unique_ptr<Impl> pimpl;
};

auto createReaction(const MineralReaction& reaction, const ChemicalSystem& system) -> Reaction;

} // namespace Reaktor
