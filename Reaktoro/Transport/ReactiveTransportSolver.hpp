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
#include <list>

// Reaktoro's includes
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumSolver.hpp>
#include <Reaktoro/Core/ChemicalOutput.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Transport/TransportSolver.hpp>

namespace Reaktoro {

struct Results{
    double smart_total;
    double smart_total_ideal_search;
    double smart_total_ideal_search_store;

    double conv_total;
};

    /// Use this class to collect modeling results per one step of reactive transport.
struct ReactiveTransportResult{

    /// Total cpu times for reactive transport and equilibrium
    double rt_time = 0.0;
    double eq_time = 0.0;

    /// Number of cells
    int ncells;
    int nsteps;

    /// Flag for the smart equilibirum
    bool is_smart;

    /// Contains reactive transport tracking results
    EquilibriumResult equilibrium;

public:

    ReactiveTransportResult(const int & ncells, const int & nsteps, const bool & is_smart);
    auto resetTime() -> void;
};

/// Use this class to postprocess and analyse accumulated results of reactive transport.
class ReactiveTransportProfiler{

    /// Name of the file and folder with a status output
    std::string folder;

    /// Used to store statistics information about the smart equilibrium algorithm.
    struct Statistics
    {
        /// Total time for search operations
        double time_estimate = 0.0;

        /// Time for search operations (part of estimation)
        double time_search = 0.0;

        /// Time for matrix-vector multiplications (part of estimation)
        double time_mat_vect_mult = 0.0;

        /// Time for acceptance test (part of estimation)
        double time_acceptance = 0.0;

        /// Total time for learn operations
        double time_learn = 0.0;

        /// Time for store operations (part of learning)
        double time_store = 0.0;

        /// Time for search operations (part of learning)
        double time_gibbs_min = 0.0;

        /// The size of the search tree
        Index tree_size = 0;

        /// Counter for the smart statuses
        int smart_counter = 0;

        /// Counter for the smart statuses
        int total_counter = 0;

    };

    /// The vector of the  CPU times for learning and estimating time
    Statistics total_stats;

    /// The vector of the  CPU times for learning and estimating time
    std::vector<Statistics> step_stats;

    /// The vector of the  CPU times for learning and estimating time
    std::vector<bool> statuses;

    /// Vector of rective transport times and equilibrium times
    std::vector<double> rt_times;
    std::vector<double> eq_times;

    /// Flag whether smart of conventional solver was used
    bool smart;

public:

    ReactiveTransportProfiler(const std::string & results_folder, const bool & smart);
    auto process(ReactiveTransportResult & rt_result) -> void;
    auto output(const std::string & file, const Index & step) -> void;
    auto console(const Index & step) -> void;
    auto summarize(Results& results) -> void;

};
/*
/// Enum type to initialize profilers
enum Profiling{
    RT = 1,
    EQ = 2,     // step-precision
    CK = 3,
    Total = 4,
    EQ_CW = 5   // cell-precision
};


/// Use this class to track, count, and document the equilibrium solver statuses
struct SmartSolverStatus{

    SmartSolverStatus(const std::string & results_folder, const std::string & file);

    /// Update the output stream with collected statuses
    auto output(const Index & i) -> void;

    /// Name of the file and folder with a status output
    std::string file, folder;

    /// The list of booleans indicating weather smart estimation was triggered
    std::vector<bool> statuses;

    /// Counter for the smart statuses
    int smart_counter = 0;

    /// Counter for the total statuses
    int total_counter = 0;

};

/// Use this class for profiling reactive transport solver components (on the step level)
struct Profiler{

    using timepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;
    using clock     = std::chrono::high_resolution_clock;

    Profiler(Profiling what);

    auto startProfiling() -> void;
    auto endProfiling() -> void;
    auto getProfilingSubject() const -> Profiling;
    auto operator==(const Profiler& p) const -> bool;

    virtual auto fileOutput(const std::string & file) -> void;

    /// Enum indicating which part of the reactive transport is profiled
    Profiling subject;

    /// The start and finish time of profiling
    timepoint start;

    /// The vector of the elapsed CPU times
    std::vector<double> times;

};

/// Use this class for profiling reactive transport components (on the cell level)
struct EquilibriumProfiler : public Profiler {

    EquilibriumProfiler(Profiling what);

    auto consoleOutput(int step) -> void;
    auto fileOutput(const std::string & file)-> void;

    /// Used to store statitics information about the smart equilibrium algorithm.
    struct Statistics
    {
        /// Total time for search operations
        double time_estimate = 0.0;

        /// Total time for learn operations
        double time_learn = 0.0;

        /// Time for search operations (part of estimation)
        double time_search = 0.0;

        /// Time for matrix-vector multiplications (part of estimation)
        double time_matrix_vector_mult = 0.0;

        /// Time for acceptance test (part of estimation)
        double time_acceptance_test = 0.0;

        /// Time for store operations (part of learning)
        double time_store = 0.0;

        /// Time for search operations (part of learning)
        double time_gibbs_minimization = 0.0;

        /// The size of the search tree
        Index tree_size = 0;
    };

    /// The vector of the  CPU times for learning and estimating time
    std::vector<Statistics> estimating_stats;
    std::vector<Statistics> learning_stats;

    /// The height of the tree storing the reference states
    int tree_height = 0;

};
*/

/// Use this class for solving reactive transport problems.
class ReactiveTransportSolver
{
public:
    /// Construct a default ReactiveTransportSolver instance.
    ReactiveTransportSolver(const ChemicalSystem& system);

    /// Set options of the reactive transport modeling
    /// TODO: make with variadic
    /// priving parameter string and value
    auto setOptions() -> void;

    /// Temporary:
    auto setSmartSolver(const bool & smart) -> void;

    auto setProfiling(const bool & profiling) -> void;

    /// Initialize the mesh discretizing the computational domain for reactive transport.
    auto setMesh(const Mesh& mesh) -> void;

    /// Initialize the velocity of the reactive transport model.
    auto setVelocity(double val) -> void;

    /// Initialize the diffusion of the reactive transport model.
    auto setDiffusionCoeff(double val) -> void;

    /// Initialize boundary conditions of the reactive transport model.
    /// Method initializes the boundary condition of the reactive transport model
    /// by setting there equilibrated ChemicalState
    /// @see ChemicalState
    auto setBoundaryState(const ChemicalState& state) -> void;

    /// Initialize time step of the reactive transport sequential algorithm.
    auto setTimeStep(double val) -> void;

    /// Get the chemical system initialized to the reactive transport model.
    /// @see ChemicalSystem
    auto system() const -> const ChemicalSystem& { return system_; }

    /// Add the output to the reactive transport modelling
    /// This method add new chemical output to the list of existing outputs
    /// and return the added element
    /// @see ChemicalOutput
    auto output() -> ChemicalOutput;

    /// Initialize components of the reactive transport model
    /// Initialize the mesh, the number of elements and cells, amounts of a fluid and
    /// solid elements as well as amounts of elements on each cell of the mesh.
    /// Finally, it initializes the matrix of the systems generated by discretized
    /// transport equation
    auto initialize() -> void;

    /// Make a step of the reactive transport time-stepping scheme
    /// @see ChemicalField
    /// auto step(ChemicalField& field) -> void;

    /// Make a step of the reactive transport time-stepping scheme
    /// @see ChemicalField
    auto step(ChemicalField& field, ReactiveTransportResult & result) -> void;

    /// Output profiling results to the file
    auto outputProfiling(const std::string & folder) -> void;

    /// Output profiling results to the console
    auto outputProfiling() -> void;

    /// Set equilibrium options to the equilibrium solvers
    /// @see EquilibriumOptions
    auto setEquilibriumOptions(const EquilibriumOptions& options) -> void;

    /// Options to customize the modelling of the reaktive transport
    struct Options{

        /// The flag indicating weather smart or conventional equilibrium solver is initialized
        bool smart = false;

        /// The flag indicating weather profiling of the reactive trasport is on or off
        bool profiling = false;
    };
    Options options;

private:
    /// The chemical system common to all degrees of freedom in the chemical field.
    ChemicalSystem system_;

    /// The solver for solving the transport equations
    TransportSolver transportsolver;

    /// The solver for solving the equilibrium equations using classical approach
    std::unique_ptr<EquilibriumSolver> equilibriumsolver;

    /// The solver for solving the equilibrium equations using smart on-demand learning algorithm
    std::unique_ptr<SmartEquilibriumSolver> smart_equilibriumsolver;

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

    /// Name of the file and folder with a status output
    std::string folder;

};

} // namespace Reaktoro