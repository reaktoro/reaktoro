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

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

/// An alias to a map of an array defining a vector of double floating-point numbers.
using VectorMap = Eigen::Map<Vector>;

/// An alias to a map of an array defining a matrix of double floating-point numbers.
using MatrixMap = Eigen::Map<Matrix>;

/// An alias to a map of an array defining a vector of indices.
using IndicesMap = Eigen::Map<Eigen::Matrix<Index, -1, 1>>;

/// A type that contains the values of a scalar field and its derivatives.
struct ChemicalField
{
    /// The values of the scalar chemical field.
    Vector values;

    /// The sensitivity of the scalar chemical field with respect to temperature.
    Vector T;

    /// The sensitivity of the scalar chemical field with respect to pressure.
    Vector P;

    /// The sensitivity of the scalar chemical field with respect to molar amounts of elements in equilibrium partition.
    Matrix be;

    /// The sensitivity of the scalar chemical field with respect to molar amounts of species in kinetic partition.
    Matrix nk;
};

/// A type that describes a solver for many chemical calculations.
class ChemicalSolver
{
public:
    /// Construct a default ChemicalSolver instance.
    ChemicalSolver();

    /// Construct a ChemicalSolver instance with given chemical system and field size.
    ChemicalSolver(const ChemicalSystem& system, unsigned size);

    /// Construct a ChemicalSolver instance with given reaction system and field size.
    ChemicalSolver(const ReactionSystem& reactions, unsigned size);

    /// Destroy a ChemicalSolver instance.
    virtual ~ChemicalSolver();

    /// Set the partitioning of the chemical system.
    auto setPartition(const Partition& partition) -> void;

    /// Set the chemical state of all points in the field.
    /// @param state The chemical state to be set in all field points.
    auto setState(const ChemicalState& state) -> void;

    /// Set the chemical state of points in the field with given indices.
    /// @param state The chemical state to be set in all selected field points.
    /// @param indices The indices of the selected field points.
    auto setState(const ChemicalState& state, const IndicesMap& indices) -> void;

    /// Return the porosity field.
    auto porosity() const -> const ChemicalField&;

    /// Return the saturation field of a fluid phase.
    /// @param ifluidphase The index of the fluid phase.
    auto saturation(unsigned ifluidphase) const -> const ChemicalField&;

    /// Return the density field of a fluid phase.
    /// @param ifluidphase The index of the fluid phase.
    auto density(unsigned ifluidphase) const -> const ChemicalField&;

    /// Equilibrate the chemical state at every field point.
    /// @param T The temperature values at every field point (in units of K).
    /// @param P The pressure values at every field point (in units of Pa).
    /// @param be The molar amounts of the equilibrium elements at every field point (in units of mol).
    auto equilibrate(const VectorMap& T, const VectorMap& P, const MatrixMap& be) -> void;

    /// React the chemical state at every field point.
    /// @param t The current time (in units of seconds)
    /// @param dt The time step to be performed (in units of seconds)
    auto react(double t, double dt) -> void;
};

} // namespace Reaktoro

