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
#include <Reaktoro/Common/ScalarTypes.hpp>
#include <Reaktoro/Thermodynamics/Reactions/MineralCatalyst.hpp>
#include <Reaktoro/Thermodynamics/Reactions/MineralMechanism.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;
class Reaction;
class ReactionEquation;

class MineralReaction
{
public:
    /// Construct a default MineralReaction instance.
    MineralReaction();

    /// Construct a MineralReaction instance with the mineral name.
    MineralReaction(std::string mineral);

    /// Set the name of the mineral species in the reaction.
    auto setMineral(std::string mineral) -> MineralReaction&;

    /// Set the equation of the mineral reaction.
    /// @see ReactionEquation
    auto setEquation(const ReactionEquation& equation) -> MineralReaction&;

    /// Set the equation of the mineral reaction.
    /// The format of @c equation can be found in ReactionEquation.
    /// @param equation The @c string defining the reaction equation
    /// @see ReactionEquation
    auto setEquation(std::string equation) -> MineralReaction&;

    /// Set the equilibrium constant of the mineral reaction (in natural log scale).
    /// If no equilibrium contant is provided, it will be calculated from the
    /// standard Gibbs energies of the species in the reaction.
    auto setEquilibriumConstant(const ThermoScalarFunction& lnk) -> MineralReaction&;

    /// Set the specific surface area of the mineral.
    /// The specific surface area of the mineral can be set using units
    /// that are convertible to either m<sup>2</sup>/g or m<sup>2</sup>/m<sup>3</sup>.
    /// Whichever unit is used, the specified specific surface area will be automatically
    /// converted to a molar surface area with units of m<sup>2</sup>/mol.
    /// @param value The value of the specific surface area
    /// @param unit The units of the specific surface area (must be convertible to either m2/g or m2/m3)
    auto setSpecificSurfaceArea(double value, std::string unit) -> MineralReaction&;

    /// Adds a mineral mechanism to the kinetic rate model of the mineral reaction
    /// @see MineralMechanism
    auto addMechanism(std::string mechanism) -> MineralReaction&;

    /// Adds a mineral mechanism to the kinetic rate model of the mineral reaction
    /// @param mechanism The mechanism to be considered in the kinetic rate model
    /// @see MineralMechanism
    auto addMechanism(const MineralMechanism& mechanism) -> MineralReaction&;

    /// Set the mineral mechanisms of the kinetic rate model of the mineral reaction.
    /// @param mechanisms The mechanisms of the kinetic rate model
    /// @see MineralMechanism
    auto setMechanisms(const std::vector<MineralMechanism>& mechanisms) -> MineralReaction&;

    /// Return the name of the mineral species in the reaction.
    auto mineral() const -> std::string;

    /// Return the equation of the mineral reaction as a @c string.
    /// @see ReactionEquation
    auto equation() const -> const ReactionEquation&;

    /// Return the equilibrium constant of the mineral reaction.
    auto equilibriumConstant() const -> const ThermoScalarFunction&;

    /// Return the specific surface area of the mineral (in units of m2/kg).
    auto specificSurfaceArea() const -> double;

    /// Return the volumetric surface area of the mineral (in units of m2/m3).
    auto volumetricSurfaceArea() const -> double;

    /// Return the mineral mechanisms of the kinetic rate model of the mineral reaction.
    /// @see MineralMechanism
    auto mechanisms() const -> const std::vector<MineralMechanism>&;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

auto createReaction(const MineralReaction& reaction, const ChemicalSystem& system) -> Reaction;

} // namespace Reaktoro
