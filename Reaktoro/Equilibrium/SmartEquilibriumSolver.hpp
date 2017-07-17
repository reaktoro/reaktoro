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

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalState;
class ChemicalSystem;
class Partition;
struct EquilibriumOptions;
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
    auto learn(ChemicalState& state, double T, double P, const Vector& be) -> EquilibriumResult;

    /// Estimate the equilibrium state using sensitivity derivatives.
    auto estimate(ChemicalState& state, double T, double P, const Vector& be) -> EquilibriumResult;

    /// Solve an equilibrium calculation either by
    auto solve(ChemicalState& state, double T, double P, const Vector& be) -> EquilibriumResult;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
