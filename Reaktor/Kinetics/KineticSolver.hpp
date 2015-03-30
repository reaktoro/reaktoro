// Reaktor is a C++ library for computational reaction modelling.
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

namespace Reaktor {

// Reaktor forward declarations
class ChemicalState;
class ReactionSystem;
class KineticProblem;
struct KineticOptions;

/// A class that represents a solver for chemical kinetics problems.
/// @see KineticProblem
class KineticSolver
{
public:
    /// Construct a default KineticSolver instance.
    KineticSolver();

    /// Construct a copy of a KineticSolver instance.
    KineticSolver(const KineticSolver& other);

    /// Destroy the KineticSolver instance.
    virtual ~KineticSolver();

    /// Assign a KineticSolver instance to this instance.
    auto operator=(KineticSolver other) -> KineticSolver&;

    /// Set the options for the chemical kinetics calculation.
    auto setOptions(const KineticOptions& options) -> void;

    /// Set the chemical kinetics problem.
    auto setProblem(const KineticProblem& problem) -> void;

    /// Integrate one step of the chemical kinetics problem.
    /// @param state The chemical state of the system
    /// @param[in,out] t The current time of the integration, updated after the calculation (in units of seconds)
    auto step(ChemicalState& state, double& t) -> void;

    /// Integrate one step of the chemical kinetics problem with a time step that does not go beyond a specified one.
    /// @param state The chemical state of the system
    /// @param[in,out] t The current time of the integration, updated after the calculation (in units of seconds)
    /// @param tfinal The final time of the integration (in units of seconds)
    auto step(ChemicalState& state, double& t, double tfinal) -> void;

    /// Solve the chemical kinetics problem from a given initial time to a final time.
    /// @param state The chemical state of the system
    /// @param t The start time of the integration (in units of seconds)
    /// @param dt The step to be used for the integration from `t` to `t + dt` (in units of seconds)
    auto solve(ChemicalState& state, double t, double dt) -> void;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktor
