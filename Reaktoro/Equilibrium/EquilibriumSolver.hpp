// Reaktoro is a C++ library for computational reaction modelling.
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
#include <memory>

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalState;
class ChemicalSystem;
class Partition;
struct EquilibriumOptions;
struct EquilibriumResult;

class EquilibriumSolver
{
public:
    /// Construct an EquilibriumSolver instance
    explicit EquilibriumSolver(const ChemicalSystem& system);

    /// Construct a copy of an EquilibriumSolver instance
    EquilibriumSolver(const EquilibriumSolver& other);

    /// Destroy this EquilibriumSolver instance
    virtual ~EquilibriumSolver();

    /// Assign a copy of an EquilibriumSolver instance
    auto operator=(EquilibriumSolver other) -> EquilibriumSolver&;

    /// Set the options of the equilibrium solver
    auto setOptions(const EquilibriumOptions& options) -> void;

    /// Set the partition of the chemical system
    auto setPartition(const Partition& partition) -> void;

    /// Set the partition of the chemical system as a formatted string
    auto setPartition(std::string partition) -> void;

    /// Find an initial feasible guess for an equilibrium problem
    /// @param state[in,out] The initial guess and the final state of the equilibrium approximation
    /// @param be The molar amounts of the elements in the equilibrium partition
    auto approximate(ChemicalState& state, const Vector& be) -> EquilibriumResult;

    /// Solve an equilibrium problem
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    /// @param be The molar amounts of the elements in the equilibrium partition
    auto solve(ChemicalState& state, const Vector& be) -> EquilibriumResult;

    /// Compute the partial derivatives @f$\left.\frac{\partial n}{\partial T}\right|_{P,b}@f$.
    /// These derivatives provide a measure of how much the equilibrium molar amounts of the species
    /// in the equilibrium partition change with an infinitesimal change in temperature. They are useful
    /// when solving non-linear problems that involve equilibrium calculations and derivatives with
    /// respect to temperature.
    auto dndT() -> Vector;

    /// Compute the partial derivatives @f$\left.\frac{\partial n}{\partial P}\right|_{T,b}@f$.
    /// These derivatives provide a measure of how much the equilibrium molar amounts of the species
    /// in the equilibrium partition change with an infinitesimal change in pressure. They are useful
    /// when solving non-linear problems that involve equilibrium calculations and derivatives with
    /// respect to pressure.
    auto dndP() -> Vector;

    /// Compute the partial derivatives @f$\left.\frac{\partial n}{\partial b}\right|_{T,P}@f$.
    /// These derivatives provide a measure of how much the equilibrium molar amounts of the species
    /// in the equilibrium partition change with an infinitesimal change in molar amounts of elements.
    /// They are useful when solving non-linear problems that involve equilibrium calculations and
    /// derivatives with respect to element molar amounts.
    auto dndb() -> Matrix;

    /// Compute the partial derivatives @f$\left.\frac{\partial n}{\partial t}\right|_{T,P}@f$.
    /// These derivatives provide a measure of how much the equilibrium molar amounts of the species
    /// in the equilibrium partition change with an infinitesimal change in time. They are useful when
    /// solving non-linear problems that involve equilibrium calculations and derivatives with respect
    /// to time.
    auto dndt(const Vector& dbdt) -> Vector;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
