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
#include <string>

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Profiling.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumSolver.hpp>
#include <Reaktoro/Kinetics/SmartKineticOptions.hpp>
#include <Reaktoro/Kinetics/SmartKineticResult.hpp>
#include <Reaktoro/Optimization/Canonicalizer.hpp>

namespace Reaktoro {

/// A class that represents a solver for chemical kinetics problems.
/// @see KineticProblem
class SmartKineticSolverBase
{
public:
    /// Construct a SmartKineticSolverBase instance.
    explicit SmartKineticSolverBase(const ReactionSystem& reactions, const Partition& partition);

    /// Destroy the SmartKineticSolverBase instance.
    virtual ~SmartKineticSolverBase();

    /// Set the partition of the chemical system.
    /// Use this method to specify the equilibrium, kinetic, and inert species.
    [[deprecated("SmartKineticSolverBase::setPartition is deprecated. Use constructor SmartKineticSolverBase(const ReactionSystem&, const Partition&) instead.")]]
    auto setPartition(const Partition& partition) -> void;

    /// Set the options for the chemical kinetics calculation.
    /// @see SmartKineticOptions
    auto setOptions(const SmartKineticOptions& options) -> void;

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

    /// Initialize the chemical kinetics solver before integration.
    /// This method should be invoked whenever the user intends to make a call to `KineticsSolver::step`.
    /// @param[in,out] state The state of the chemical system
    /// @param tstart The start time of the integration.
    auto initialize(ChemicalState& state, double tstart) -> void;

    /// Initialize the chemical kinetics solver before integration
    /// with the provided vector of unknowns `benk` = [be, nk].
    /// This method should be invoked whenever the user intends to make a call to `KineticsSolver::solve`.
    /// @param[in,out] state The state of the chemical system
    /// @param tstart The start time of the integration
    /// @param benk The initial vector of unknowns `benk` = [be, nk].
    auto initialize(ChemicalState& state, double tstart, VectorConstRef benk) -> void;

    /// Function representing the right-hand side of the system od ODEs to be evaluated.
    /// @param[in,out] state The state of the chemical system that provides chemical properties needed for the evaluation
    /// @param t The current time of the right-hand side evaluation
    /// @param res The result of the evaluation
    auto function(ChemicalState& state, double t, VectorConstRef u, VectorRef res) -> int;

    /// Function representing the Jacobian of the right-hand side of the system od ODEs to be evaluated.
    /// @param[in,out] state The state of the chemical system that provides chemical properties needed for the evaluation
    /// @param t The current time of the right-hand side evaluation
    /// @param res The result of the evaluation
    auto jacobian(ChemicalState& state, double t, VectorConstRef u, MatrixRef res) -> int;

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

    /// Learn how to perform a full equilibrium calculation (with tracking)
    virtual auto learn(ChemicalState& state, double& t, double dt) -> void = 0;

    /// Estimate the equilibrium state using sensitivity derivatives (profiling the expenses)
    virtual auto estimate(ChemicalState& state, double& t, double dt) -> void = 0;

    /// Return the result of the last smart kinetic calculation.
    /// @see SmartKineticResult
    auto result() const -> const SmartKineticResult&;

    // Return properties of the chemical state provided by the KineticSolver.
    /// @see ChemicalProperties
    auto properties() const -> const ChemicalProperties&;

    // Return reaction system of the chemical state provided by the KineticSolver.
    /// @see ReactionSystem
    auto reactions() const -> const ReactionSystem&;

    // Return partition of the chemical state provided by the KineticSolver.
    /// @see Partition
    auto partition() const -> const Partition&;

    /// Output clusters created during the ODML algorithm.
    virtual auto outputInfo() const -> void = 0;

protected:

    /// The kinetically-controlled chemical reactions
    ReactionSystem _reactions;

    /// The chemical system instance
    ChemicalSystem system;

    /// The partition of the species in the chemical system
    Partition _partition;

    /// The options of the kinetic solver
    SmartKineticOptions options;

    /// The result of the smart kinetic solver
    SmartKineticResult _result;

    // The solver for solving the equilibrium equations using classical approach
    EquilibriumSolver equilibrium;

    /// The solver for solving the equilibrium equations using smart on-demand learning algorithm
    SmartEquilibriumSolver smart_equilibrium;

    /// The sensitivity of the equilibrium state
    EquilibriumSensitivity sensitivity;

    /// The chemical properties of the chemical system
    ChemicalProperties _properties;

    /// The ODE solver instance
    ODESolver ode;

    /// The canonicalizer used to determine primary and secondary species
    Canonicalizer canonicalizer;

    /// The indices of the equilibrium, kinetic, and inert species
    Indices ies, iks, iis;

    /// The indices of the elements in the equilibrium and kinetic partition
    Indices iee, ike;

    /// The number of equilibrium, kinetic, inert, and all species
    Index Ne, Nk, Ni, N;

    /// The number of elements in the equilibrium partition
    Index Ee, E;

    /// The formula matrix of the equilibrium species
    Matrix Ae, Ak;

    /// The stoichiometric matrix w.r.t. the equilibrium and kinetic species
    Matrix Se, Sk;

    /// The coefficient matrix `A` of the kinetic rates
    Matrix A;

    /// The coefficient matrix `B` of the source rates
    Matrix B;

    /// The temperature of the chemical system (in units of K)
    double T;

    /// The pressure of the chemical system (in units of Pa)
    double P;

    /// The amounts of the species in the chemical system
    Vector n;

    /// The molar composition of the equilibrium and kinetic species
    Vector ne, nk;

    /// The molar abundance of the elements in the equilibrium species
    Vector be;

    /// The combined vector of elemental molar abundance and composition of kinetic species [be nk]
    Vector benk;

    /// The vector with the values of the reaction rates
    ChemicalVector rates;

    /// The partial derivatives of the reaction rates `rates` w.r.t. to `be`, `ne`, `nk`, `and `u = [be nk]`
    Matrix drdbe, drdne, drdnk, drdu;

    /// The source term
    ChemicalVector q;

    /// The partial derivatives of the source rates `q` w.r.t. to `be`, `ne`, `nk`, `and `u = [be nk]`
    Matrix dqdbe, dqdne, dqdnk, dqdu;

    /// The function that calculates the source term in the problem
    std::function<ChemicalVector(const ChemicalProperties&)> source_fn;

    /// The sensitivity matrix of combined vector of elemental molar abundance and composition of kinetic species [be nk]
    Matrix benk_S;

};

} // namespace Reaktoro
