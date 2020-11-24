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
#include <string>

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

// Forward declarations (classes)
class ChemicalState;
class ChemicalProperties;
class Partition;
class ReactionSystem;

// Forward declarations (structs)
struct KineticOptions;
struct KineticResult;

/// A class that represents a solver for chemical kinetics problems.
/// @see KineticProblem
class KineticSolver
{
public:
    /// Construct a default KineticSolver instance.
    [[deprecated("KineticSolver() is deprecated. Use constructor KineticSolver(const ReactionSystem&, const Partition&) instead.")]]
    KineticSolver();

    /// Construct a KineticSolver instance
    [[deprecated("KineticSolver(const ReactionSystem&) is deprecated. Use constructor KineticSolver(const ReactionSystem&, const Partition&) instead.")]]
    explicit KineticSolver(const ReactionSystem& reactions);

    /// Construct a KineticSolver instance.
    KineticSolver(const ReactionSystem& reactions, const Partition& partition);

    /// Assign a KineticSolver instance to this instance.
    auto operator=(KineticSolver other) -> KineticSolver&;

    /// Destroy the KineticSolver instance.
    virtual ~KineticSolver();

    /// Set the options for the chemical kinetics calculation.
    auto setOptions(const KineticOptions& options) -> void;

    /// Set the partition of the chemical system.
    /// Use this method to specify the equilibrium, kinetic, and inert species.
    [[deprecated("KineticSolver::setPartition is deprecated. Use constructor KineticSolver(const ReactionSystem&, const Partition&) instead.")]]
    auto setPartition(const Partition& partition) -> void;

    /// Add a source to the chemical kinetics problem.
    /// @param state The chemical state representing the source.
    /// @param volumerate The volumetric rate of the source.
    /// @param units The units of the volumetric rate (compatible with m3/s).
    auto addSource(const ChemicalState& state, double volumerate, const std::string& units) -> void;

    /// Add a phase sink to the chemical kinetics problem.
    /// This method allows the chemical kinetics problem to account for
    /// the sink (i.e., the removal) of a phase from the system.
    /// @param phase The name of the phase.
    /// @param volumerate The volumetric rate of the phase removal.
    /// @param units The units of the volumetric rate (compatible with m3/s).
    auto addPhaseSink(const std::string& phase, double volumerate, const std::string& units) -> void;

    /// Add a fluid sink to the chemical kinetics problem.
    /// This method allows the chemical kinetics problem to account for
    /// the sink (i.e., the removal) of fluid from the system.
    /// @param volumerate The volumetric rate of the fluid removal.
    /// @param units The units of the volumetric rate (compatible with m3/s).
    auto addFluidSink(double volumerate, const std::string& units) -> void;

    /// Add a solid sink to the chemical kinetics problem.
    /// This method allows the chemical kinetics problem to account for
    /// the sink (i.e., the removal) of solid from the system.
    /// @param volumerate The volumetric rate of the solid removal.
    /// @param units The units of the volumetric rate (compatible with m3/s).
    auto addSolidSink(double volumerate, const std::string& units) -> void;

    /// Initialize the chemical kinetics solver before integration.
    /// This method should be invoked whenever the user intends to make a call to `KineticsSolver::step`.
    /// @param state The state of the chemical system
    /// @param tstart The start time of the integration.
    auto initialize(ChemicalState& state, double tstart) -> void;

    /// Initialize the chemical kinetics solver before integration
    /// with the provided vector of unknowns `benk` = [be, nk].
    /// This method should be invoked whenever the user intends to make a call to `KineticsSolver::solve`.
    /// @param state The state of the chemical system
    /// @param tstart The start time of the integration
    /// @param benk The initial vector of unknowns `benk` = [be, nk].
    auto initialize(ChemicalState& state, double tstart, VectorConstRef benk) -> void;

    /// Integrate one step of the chemical kinetics problem.
    /// @param state The kinetic state of the system
    /// @param[in,out] t The current time of the integration (in units of seconds)
    /// @return The updated current time after the kinetic step.
    auto step(ChemicalState& state, double t) -> double;

    /// Integrate one step of the chemical kinetics problem with a time step 'dt' that does not go beyond a specified
    /// time `tfinal`, i.e., `t + dt` < `tfinal`. Here, the time step `dt` is determined by the CVODE solver.
    /// @param state The kinetic state of the system
    /// @param[in,out] t The current time of the integration (in units of seconds)
    /// @param tfinal The final time of the integration (in units of seconds)
    /// @return The updated current time after the kinetic step.
    auto step(ChemicalState& state, double t, double tfinal) -> double;

    /// Solve the chemical kinetics problem from a given initial time to a final time.
    /// Used in function of KineticPath class, where the uniform stepping in time is required.
    /// @param state The kinetic state of the system
    /// @param t The start time of the integration (in units of seconds)
    /// @param dt The step to be used for the integration from `t` to `t + dt` (in units of seconds)
    auto solve(ChemicalState& state, double t, double dt) -> double;

    /// Solve the chemical kinetics problem from a given initial time to a final time with provided amount of elements
    /// updated from the transport.
    /// Used in reactive transport solver, which provides updated element amounts.
    /// @param state The kinetic state of the system
    /// @param t The start time of the integration (in units of seconds)
    /// @param dt The step to be used for the integration from `t` to `t + dt` (in units of seconds)
    /// @param b The amount of elements updated from the transport
    auto solve(ChemicalState& state, double t, double dt, VectorConstRef b) -> void;

    /// Return the result of the last kinetic calculation.
    auto result() const -> const KineticResult&;

    /// Return properties of the chemical state provided by the KineticSolver
    auto properties() const -> const ChemicalProperties&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
