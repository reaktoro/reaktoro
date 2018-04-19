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

namespace Reaktoro {

// Forward declarations
class ChemicalState;
class ChemicalState;
class Partition;
class ReactionSystem;
struct KineticOptions;

/// A class that represents a solver for chemical kinetics problems.
/// @see KineticProblem
class KineticSolver
{
public:
    /// Construct a default KineticSolver instance.
    KineticSolver();

    /// Construct a KineticSolver instance.
    explicit KineticSolver(const ReactionSystem& reactions);

    /// Construct a copy of a KineticSolver instance.
    KineticSolver(const KineticSolver& other);

    /// Destroy the KineticSolver instance.
    virtual ~KineticSolver();

    /// Assign a KineticSolver instance to this instance.
    auto operator=(KineticSolver other) -> KineticSolver&;

    /// Set the options for the chemical kinetics calculation.
    auto setOptions(const KineticOptions& options) -> void;

    /// Set the partition of the chemical system.
    /// Use this method to specify the equilibrium, kinetic, and inert species.
    auto setPartition(const Partition& partition) -> void;

    /// Add a source to the chemical kinetics problem.
    /// @param state The chemical state representing the source.
    /// @param volumerate The volumetric rate of the source.
    /// @param units The units of the volumetric rate (compatible with m3/s).
    auto addSource(const ChemicalState& state, double volumerate, std::string units) -> void;

    /// Add a phase sink to the chemical kinetics problem.
    /// This method allows the chemical kinetics problem to account for
    /// the sink (i.e., the removal) of a phase from the system.
    /// @param phase The name of the phase.
    /// @param volumerate The volumetric rate of the phase removal.
    /// @param units The units of the volumetric rate (compatible with m3/s).
    auto addPhaseSink(std::string phase, double volumerate, std::string units) -> void;

    /// Add a fluid sink to the chemical kinetics problem.
    /// This method allows the chemical kinetics problem to account for
    /// the sink (i.e., the removal) of fluid from the system.
    /// @param volumerate The volumetric rate of the fluid removal.
    /// @param units The units of the volumetric rate (compatible with m3/s).
    auto addFluidSink(double volumerate, std::string units) -> void;

    /// Add a solid sink to the chemical kinetics problem.
    /// This method allows the chemical kinetics problem to account for
    /// the sink (i.e., the removal) of solid from the system.
    /// @param volumerate The volumetric rate of the solid removal.
    /// @param units The units of the volumetric rate (compatible with m3/s).
    auto addSolidSink(double volumerate, std::string units) -> void;

    /// Initialize the chemical kinetics solver before .
    /// This method should be invoked whenever the user intends to make a call to `KineticsSolver::step`.
    /// @param state The state of the chemical system
    /// @param tstart The start time of the integration.
    auto initialize(ChemicalState& state, double tstart) -> void;

    /// Integrate one step of the chemical kinetics problem.
    /// @param state The kinetic state of the system
    /// @param[in,out] t The current time of the integration (in units of seconds)
    /// @return The updated current time after the kinetic step.
    auto step(ChemicalState& state, double t) -> double;

    /// Integrate one step of the chemical kinetics problem with a time step that does not go beyond a specified one.
    /// @param state The kinetic state of the system
    /// @param[in,out] t The current time of the integration (in units of seconds)
    /// @param tfinal The final time of the integration (in units of seconds)
    /// @return The updated current time after the kinetic step.
    auto step(ChemicalState& state, double t, double tfinal) -> double;

    /// Solve the chemical kinetics problem from a given initial time to a final time.
    /// @param state The kinetic state of the system
    /// @param t The start time of the integration (in units of seconds)
    /// @param dt The step to be used for the integration from `t` to `t + dt` (in units of seconds)
    auto solve(ChemicalState& state, double t, double dt) -> void;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
