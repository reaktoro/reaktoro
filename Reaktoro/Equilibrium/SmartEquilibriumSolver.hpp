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

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProperties;
class ChemicalState;
class ChemicalSystem;
class Partition;
struct EquilibriumOptions;
class EquilibriumProblem;
struct EquilibriumResult;

/// A class used to perform equilibrium calculations using machine learning scheme.
class SmartEquilibriumSolver
{
public:
    /// Construct a default SmartEquilibriumSolver instance.
    SmartEquilibriumSolver();

    /// Construct an SmartEquilibriumSolver instance
    explicit SmartEquilibriumSolver(const ChemicalSystem& system);

    /// Construct a copy of an SmartEquilibriumSolver instance.
    SmartEquilibriumSolver(const SmartEquilibriumSolver& other);

    /// Assign an SmartEquilibriumSolver instance to this.
    auto operator=(SmartEquilibriumSolver other) -> SmartEquilibriumSolver&;

    /// Destroy this SmartEquilibriumSolver instance.
    virtual ~SmartEquilibriumSolver();

    /// Set the options for the equilibrium calculation.
    auto setOptions(const EquilibriumOptions& options) -> void;

    /// Set the partition of the chemical system.
    auto setPartition(const Partition& partition) -> void;

    /// Learn how to perform a full equilibrium calculation.
    auto learn(ChemicalState& state, double T, double P, VectorXrConstRef be) -> EquilibriumResult;

    /// Learn how to perform a full equilibrium calculation.
    auto learn(ChemicalState& state, const EquilibriumProblem& problem) -> EquilibriumResult;

    /// Estimate the equilibrium state using sensitivity derivatives.
    auto estimate(ChemicalState& state, double T, double P, VectorXrConstRef be) -> EquilibriumResult;

    /// Estimate the equilibrium state using sensitivity derivatives.
    auto estimate(ChemicalState& state, const EquilibriumProblem& problem) -> EquilibriumResult;

    /// Solve an equilibrium calculation either by
    auto solve(ChemicalState& state, double T, double P, VectorXrConstRef be) -> EquilibriumResult;

    /// Solve an equilibrium problem with given equilibrium problem.
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    /// @param problem The equilibrium problem with given temperature, pressure, and element amounts.
    auto solve(ChemicalState& state, const EquilibriumProblem& problem) -> EquilibriumResult;

    /// Return the chemical properties of the calculated equilibrium state.
    auto properties() const -> const ChemicalProperties&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
