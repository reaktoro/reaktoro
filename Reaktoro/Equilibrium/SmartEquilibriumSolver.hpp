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

// C++ includes
#include <memory>

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

// Forward declarations (classes)
class ChemicalProperties;
class ChemicalState;
class ChemicalSystem;
class EquilibriumProblem;
class EquilibriumSensitivity;
class Partition;

// Forward declarations (structs)
struct SmartEquilibriumOptions;
struct SmartEquilibriumProfiling;
struct SmartEquilibriumResult;
struct SmartEquilibriumResultDuringEstimate;
struct SmartEquilibriumResultDuringLearning;

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
    auto setOptions(const SmartEquilibriumOptions& options) -> void;

    /// Set the partition of the chemical system.
    auto setPartition(const Partition& partition) -> void;

    /// Learn how to perform a full equilibrium calculation.
    auto learn(ChemicalState& state, double T, double P, VectorConstRef be) -> SmartEquilibriumResultDuringLearning;

    /// Learn how to perform a full equilibrium calculation.
    auto learn(ChemicalState& state, const EquilibriumProblem& problem) -> SmartEquilibriumResultDuringLearning;

    /// Estimate the equilibrium state using sensitivity derivatives.
    auto estimate(ChemicalState& state, double T, double P, VectorConstRef be) -> SmartEquilibriumResultDuringEstimate;

    /// Estimate the equilibrium state using sensitivity derivatives.
    auto estimate(ChemicalState& state, const EquilibriumProblem& problem) -> SmartEquilibriumResultDuringEstimate;

    /// Solve an equilibrium state.
    /// Solving of the SmartEquilibriumSolver consists of two possible stages:
    /// learning (i.e., triggering the convention approach of Gibbs' minimization problem)
    /// or estimating (i.e., smart prediction of the new states using sensitivity derivatives)
    /// of the chemical state
    auto solve(ChemicalState& state, double T, double P, VectorConstRef be) -> SmartEquilibriumResult;

    /// Solve an equilibrium problem with given equilibrium problem.
    /// @param state[in,out] The initial guess and the final state of the equilibrium calculation
    /// @param problem The equilibrium problem with given temperature, pressure, and element amounts.
    auto solve(ChemicalState& state, const EquilibriumProblem& problem) -> SmartEquilibriumResult;

    /// Return the chemical properties of the calculated equilibrium state.
    auto properties() const -> const ChemicalProperties&;

    /// Return the profiling information of the operations during a smart equilibrium calculation.
    auto profiling() const -> const SmartEquilibriumProfiling&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
