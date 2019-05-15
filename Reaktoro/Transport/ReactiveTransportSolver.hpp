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
#include <chrono>
#include <iostream>
#include <fstream>

// Reaktoro's includes
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumSolver.hpp>
#include <Reaktoro/Core/ChemicalOutput.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Transport/TransportSolver.hpp>

namespace Reaktoro {

// Use this enum type to initialize profilers
enum Profiling{
    RT = 1,
    EQ = 2,
    CK = 3,
    Total = 4,
};

struct SolverStatus{

    SolverStatus(const std::string & results_folder, const std::string & file);

    /// Update the output stream with collected statuses
    auto output(const Index & i) -> void;

    // Name of the file with a status output
    std::string file;
    // Name of the folder with a status output
    std::string folder;

    /// The list of booleans indicating weather smart estimation was triggered
    std::vector<bool> statuses;
};

/// Use this class for profiling reactive transport components
class Profiler{

public:
    Profiler(Profiling what);
    auto startProfiling() -> void;
    auto endProfiling() -> void;
    auto fileOutput(const std::string & file)-> void;
    auto consoleOutput()-> void;
    auto getProfilingSubject() const -> Profiling;
    auto operator==(const Profiler& p) const -> bool;

private:
    /// Enum indicating which part of the reactive transport is profiled
    Profiling subject;

    /// The start and finish time of profiling
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> finish;

    /// The vector of the elapsed CPU times
    std::vector<double> times;

    /// The start and finish time of profiling
    std::chrono::seconds start_;
    std::chrono::seconds finish_;

    /// The vector of the elapsed CPU times
    std::vector<std::chrono::seconds> times_;

};

/// Use this class for solving reactive transport problems.
class ReactiveTransportSolver
{
public:
    /// Construct a default ReactiveTransportSolver instance.
    ReactiveTransportSolver(const ChemicalSystem& system, const bool& is_smart = false);

    auto setMesh(const Mesh& mesh) -> void;

    auto setVelocity(double val) -> void;

    auto setDiffusionCoeff(double val) -> void;

    auto setBoundaryState(const ChemicalState& state) -> void;

    auto setTimeStep(double val) -> void;

    auto system() const -> const ChemicalSystem& { return system_; }

    auto output() -> ChemicalOutput;

    auto initialize() -> void;

    auto step(ChemicalField& field) -> void;

    auto step_tracked(ChemicalField& field) -> void;

    auto profile(Profiling what) -> Profiler;

    auto trackStatus(const std::string & folder, const std::string & file) -> SolverStatus;

    auto outputProfiling(const std::string & folder) -> void;

    auto setEquilibriumOptions(const EquilibriumOptions& options) -> void;

private:
    /// The chemical system common to all degrees of freedom in the chemical field.
    ChemicalSystem system_;

    /// The solver for solving the transport equations
    TransportSolver transportsolver;

    /// The solver for solving the equilibrium equations using classical approach
    EquilibriumSolver equilibriumsolver;

    /// The solver for solving the equilibrium equations using smart on-demand learning algorithm
    SmartEquilibriumSolver smart_equilibriumsolver;

    /// The bool parameter indicating weather smart equilibrium solver is initialized
    bool smart;

    /// The list of chemical output objects
    std::vector<ChemicalOutput> outputs;

    /// The amounts of fluid elements on the boundary.
    Vector bbc;

    /// The amounts of a fluid element on each cell of the mesh.
    Matrix bf;

    /// The amounts of a solid element on each cell of the mesh.
    Matrix bs;

    /// The amounts of an element on each cell of the mesh.
    Matrix b;

    /// The current number of steps in the solution of the reactive transport equations.
    Index steps = 0;

    /// The classes to profile reactive transport computations
    std::vector<Profiler> profilers;

    /// The smart solver tracker
    std::vector<SolverStatus> status_trackers;
};

}