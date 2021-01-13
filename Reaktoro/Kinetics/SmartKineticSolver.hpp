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
#include <Reaktoro/Equilibrium/SmartEquilibriumOptions.hpp>
#include <Reaktoro/Kinetics/SmartKineticSolverBase.hpp>

namespace Reaktoro {

/// A class used to perform kinetic calculations using machine learning scheme.
class SmartKineticSolver
{
public:
    /// Construct a SmartKineticSolver instance with given reaction system and partition.
    explicit SmartKineticSolver(const ReactionSystem& reactions, const Partition& partition);

    /// Destroy this SmartKineticSolverBase instance.
    virtual ~SmartKineticSolver();

    /// Set the options for the chemical kinetics calculation.
    /// @see SmartKineticOptions
    auto setOptions(const SmartKineticOptions& options) -> void;

    /// Initialize the chemical kinetics solver before integration.
    /// This method should be invoked whenever the user intends to make a call to `KineticsSolver::step`.
    /// @param[in,out] state The state of the chemical system
    /// @param tstart The start time of the integration.
    auto initialize(ChemicalState& state, double tstart) -> void;

    /// Add a source to the chemical kinetics problem.
    /// @param state The chemical state representing the source.
    /// @param volumerate The volumetric rate of the source.
    /// @param units The units of the volumetric rate (compatible with m3/s).
    auto addSource(ChemicalState state, double volumerate, const std::string& units) -> void;

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

    /// Solve the chemical kinetics problem from a given initial time to a final time (used in KineticPath).
    /// @param[in,out] state The kinetic state of the system
    /// @param t The start time of the integration (in units of seconds)
    /// @param dt The step to be used for the integration from `t` to `t + dt` (in units of seconds)
    /// @see KineticPath
    auto solve(ChemicalState& state, double t, double dt) -> double;

    /// Solve the chemical kinetics problem from a given initial time to a final time with the given elements
    /// amount vector (used in ReactiveTransportSolver).
    /// Used in reactive transport solver, which provides updated element amounts.
    /// @param[in,out] state The kinetic state of the system
    /// @param t The start time of the integration (in units of seconds)
    /// @param dt The step to be used for the integration from `t` to `t + dt` (in units of seconds)
    /// @param b The amount of elements updated from the transport
    /// @see ReactiveTransportSolver
    auto solve(ChemicalState& state, double t, double dt, VectorConstRef b) -> void;

    /// Return the result of the last smart kinetic calculation.
    /// @see SmartKineticResult
    auto result() const -> const SmartKineticResult&;

    /// Output into about the ODML algorithm
    auto outputInfo() const -> void;

private:
    // Pointer to the smart kinetics solver instance
    std::unique_ptr<SmartKineticSolverBase> solverptr;

    // Flag to indicate whether the solver has been initialized
    bool initialized = false;

};

} // namespace Reaktoro
