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
#include <vector>

namespace Reaktoro {

// Forward declarations
class ChemicalOutput;
class ChemicalPlot;
class ChemicalSystem;
class ChemicalState;
class ChemicalState;
class Partition;
class ReactionSystem;
struct KineticOptions;

/// A class that conveniently solves kinetic path calculations.
class KineticPath
{
public:
    /// Construct a KineticPath instance.
    [[deprecated("KineticPath(const ReactionSystem&) is deprecated. Use constructor KineticPath(const ReactionSystem&, const Partition&) instead.")]]
    explicit KineticPath(const ReactionSystem& reactions);

    /// Construct a KineticPath instance.
    KineticPath(const ReactionSystem& reactions, const Partition& partition);

    /// Construct a copy of a KineticPath instance.
    KineticPath(const KineticPath& other) = delete;

    /// Destroy the KineticPath instance.
    virtual ~KineticPath();

    /// Assign a KineticPath instance to this instance.
    auto operator=(KineticPath other) -> KineticPath&;

    /// Set the options for the chemical kinetics calculation.
    auto setOptions(const KineticOptions& options) -> void;

    /// Set the partition of the chemical system.
    /// Use this method to specify the equilibrium, kinetic, and inert species.
    [[deprecated("KineticPath::setPartition is deprecated. Use constructor KineticPath(const ReactionSystem&, const Partition&) instead.")]]
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

    /// Solve the kinetic path problem (adaptively within the interval `[t0, t1]` with CVODE).
    /// @param state The initial state of the chemical system
    /// @param t0 The initial time of the kinetic path
    /// @param t1 The final time of the kinetic path
    /// @param units The time units of `t0` and `t1` (e.g., `s`, `minute`, `day`, `year`, etc.).
    auto solve(ChemicalState& state, double t0, double t1, const std::string& units = "s") -> void;

    /// Solve the kinetic path problem (uniformly with the `n` time steps of the size `dt` starting from the time `t0`,
    /// where each time step `[t0, t0 + dt]` is solved with internal CVODE-solver adaptively).
    /// @param state The initial state of the chemical system
    /// @param t0 The initial time of the kinetic path
    /// @param t1 The final time of the kinetic path
    /// @param units The time units of `t0` and `t1` (e.g., `s`, `minute`, `day`, `year`, etc.).
    auto solve(ChemicalState& state, double t0, double dt, int n, const std::string& units) -> void;

    /// Return a ChemicalPlot instance.
    /// The returned ChemicalOutput instance must be properly configured
    /// before the method EquilibriumPath::solve is called.
    /// Changes in this ChemicalOutput instance are observed by the
    /// EquilibriumPath object.
    auto output() -> ChemicalOutput;

    /// Return a ChemicalPlot instance.
    /// The returned ChemicalPlot instance must be properly configured
    /// before the method EquilibriumPath::solve is called.
    /// Changes in this ChemicalPlot instance are observed by the
    /// EquilibriumPath object.
    auto plot() -> ChemicalPlot;

    /// Return a collection of ChemicalPlot instances.
    /// The returned ChemicalPlot instances must be properly configured
    /// before the method EquilibriumPath::solve is called.
    /// Changes in theses ChemicalPlot instances are observed by the
    /// EquilibriumPath object.
    auto plots(unsigned num) -> std::vector<ChemicalPlot>;

    /// Return the chemical system in the kinetic path definition.
    auto system() const -> const ChemicalSystem&;

    /// Return the reactions in the kinetic path definition.
    auto reactions() const -> const ReactionSystem&;

    /// Return the partition of the chemical system in the kinetic path definition.
    auto partition() const -> const Partition&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
