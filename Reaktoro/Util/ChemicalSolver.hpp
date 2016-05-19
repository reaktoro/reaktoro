// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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
#include <memory>
#include <string>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalField;
class ChemicalSystem;
class ChemicalState;
class ReactionSystem;
class Partition;

/// A type that describes a solver for many chemical calculations.
class ChemicalSolver
{
public:
    /// Construct a default ChemicalSolver instance.
    ChemicalSolver();

    /// Construct a ChemicalSolver instance with given chemical system and field number of points.
    ChemicalSolver(const ChemicalSystem& system, Index npoints);

    /// Construct a ChemicalSolver instance with given reaction system and field number of points.
    ChemicalSolver(const ReactionSystem& reactions, Index npoints);

    /// Set the partitioning of the chemical system.
    auto setPartition(const Partition& partition) -> void;

    /// Set the chemical state of all points in the field.
    /// @param state The chemical state to be set in all field points.
    auto setState(const ChemicalState& state) -> void;

    /// Set the chemical state of points in the field with given indices.
    /// @param state The chemical state to be set in all selected field points.
    /// @param indices The indices of the selected field points.
    auto setStates(const ChemicalState& state, const Indices& indices) -> void;

    /// Return the chemical state at given index.
    auto state(Index i) const -> const ChemicalState&;

    /// Return the chemical states at all field points.
    auto states() const -> const std::vector<ChemicalState>&;

    /// Equilibrate the chemical state at every field point.
    /// @param T The temperature values at every field point (in units of K).
    /// @param P The pressure values at every field point (in units of Pa).
    /// @param be The molar amounts of the equilibrium elements at every field point (in units of mol).
    auto equilibrate(const double* T, const double* P, const double* be) -> void;

    /// React the chemical state at every field point.
    /// @param t The current time (in units of seconds)
    /// @param dt The time step to be performed (in units of seconds)
    auto react(double t, double dt) -> void;

    /// Return the molar amounts of the equilibrium species and their derivatives at every field point.
    auto ne() -> const std::vector<ChemicalField>&;

    /// Return the porosity at every field point.
    auto porosity() -> const ChemicalField&;

    /// Return the saturation of each fluid phase at every field point.
    auto saturations() -> const std::vector<ChemicalField>&;

    /// Return the density of each fluid phase at every field point.
    auto densities() -> const std::vector<ChemicalField>&;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro

