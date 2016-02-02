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
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalState;
class ChemicalSystem;
class EquilibriumInverseProblem;
class Partition;
struct EquilibriumOptions;
struct EquilibriumResult;

/// A class used for solving inverse equilibrium problems.
class EquilibriumInverseSolver
{
public:
    /// Construct an EquilibriumInverseSolver instance
    explicit EquilibriumInverseSolver(const ChemicalSystem& system);

    /// Construct a copy of an EquilibriumInverseSolver instance
    EquilibriumInverseSolver(const EquilibriumInverseSolver& other);

    /// Destroy this EquilibriumInverseSolver instance
    virtual ~EquilibriumInverseSolver();

    /// Assign a copy of an EquilibriumInverseSolver instance
    auto operator=(EquilibriumInverseSolver other) -> EquilibriumInverseSolver&;

    /// Set the options of the equilibrium solver
    auto setOptions(const EquilibriumOptions& options) -> void;

    /// Set the partition of the chemical system
    auto setPartition(const Partition& partition) -> void;

    /// Set the partition of the chemical system as a formatted string
    auto setPartition(std::string partition) -> void;

    /// Solve an inverse equilibrium problem.
    /// @param state[in,out] The initial guess and the final state of the inverse equilibrium calculation.
    /// @param problem The inverse equilibrium problem.
    auto solve(ChemicalState& state, const EquilibriumInverseProblem& problem) -> EquilibriumResult;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
