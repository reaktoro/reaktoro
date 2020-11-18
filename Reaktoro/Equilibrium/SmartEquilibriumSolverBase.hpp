// Reaktoro is a unified framework for modeling chemically reactive systems.
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

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Profiling.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>
#include <Reaktoro/Optimization/Canonicalizer.hpp>

namespace Reaktoro {

/// A class used to perform equilibrium calculations using machine learning scheme.
class SmartEquilibriumSolverBase
{
public:
    /// Construct an SmartEquilibriumSolverBase instance with given chemical system
    explicit SmartEquilibriumSolverBase(const ChemicalSystem& system);

    /// Construct an SmartEquilibriumSolverBase instance with given partition
    explicit SmartEquilibriumSolverBase(const Partition& partition);

    /// Destroy this SmartEquilibriumSolverBase instance.
    virtual ~SmartEquilibriumSolverBase();

    /// Set the partition of the chemical system
    [[deprecated("EquilibriumSolver::setPartition is deprecated. Use constructor EquilibriumSolver(const Partition&) instead.")]]
    auto setPartition(const Partition& partition) -> void;

    /// Set the options for the equilibrium calculation.
    auto setOptions(const SmartEquilibriumOptions& options) -> void;

    /// Solve a chemical equilibrium problem.
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    /// @param be The amounts of the elements in the equilibrium partition
    auto solve(ChemicalState& state, double T, double P, VectorConstRef be) -> SmartEquilibriumResult;

    /// Solve a chemical equilibrium problem.
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    /// @param problem The equilibrium problem with given temperature, pressure, and element amounts.
    auto solve(ChemicalState& state, const EquilibriumProblem& problem) -> SmartEquilibriumResult;

    /// Learn how to perform a full equilibrium calculation (with tracking)
    virtual auto learn(ChemicalState& state, double T, double P, VectorConstRef be) -> void = 0;

    /// Estimate the equilibrium state using sensitivity derivatives (profiling the expences)
    virtual auto estimate(ChemicalState& state, double T, double P, VectorConstRef be) -> void = 0;

    /// Return the chemical properties of the calculated equilibrium state.
    auto getProperties() const -> const ChemicalProperties&;

    /// Return the result of the last smart equilibrium calculation.
    auto getResult() const -> const SmartEquilibriumResult&;

protected:
    /// The chemical system instance
    ChemicalSystem system;

    /// The partition of the chemical system
    Partition partition;

    /// The canonicalizer used to determine primary and secondary species
    Canonicalizer canonicalizer;

    /// The options for the smart equilibrium calculation
    SmartEquilibriumOptions options;

    /// The result of the last smart equilibrium calculation.
    SmartEquilibriumResult result;

    /// The solver for the equilibrium calculations
    EquilibriumSolver solver;

    /// The chemical properties of the chemical system
    ChemicalProperties properties;

    /// The amounts of the species in the chemical system
    Vector n;

    /// The amounts of the equilibrium species
    Vector ne;

    /// The amounts of the elements in the equilibrium partition
    Vector be;

    /// The indices of the equilibrium species
    Indices ies;

    /// The indices of the equilibrium elements
    Indices iee;
};

} // namespace Reaktoro
