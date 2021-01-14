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

#include "KineticPath.hpp"

// C++ includes
#include <ostream>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalOutput.hpp>
#include <Reaktoro/Core/ChemicalPlot.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Kinetics/KineticSolver.hpp>
#include <Reaktoro/Kinetics/KineticPathProfiler.hpp>
#include <Reaktoro/Kinetics/KineticPathOptions.hpp>
#include <Reaktoro/Kinetics/SmartKineticSolver.hpp>
#include <Reaktoro/Common/Json.hpp>
#include <Reaktoro/Serialization/JsonKineticPath.hpp>

namespace Reaktoro {

struct KineticPath::Impl
{
    /// The kinetically-controlled chemical reactions
    ReactionSystem reactions;

    /// The chemical system instance
    ChemicalSystem system;

    /// The partition of the species in the chemical system
    Partition partition;

    /// The kinetic solver instance
    KineticSolver kinetic_solver;

    /// The smart kinetic solver instance
    SmartKineticSolver smart_kinetic_solver;

    /// The options of the kinetic path
    KineticPathOptions options;

    /// The output instance of the kinetic path calculation
    ChemicalOutput output;

    /// The plots of the kinetic path calculation
    std::vector<ChemicalPlot> plots;

    /// Construct an instance of KineticPath::Impl
    Impl(const ReactionSystem& reactions, const Partition& partition)
    : reactions(reactions),
      system(partition.system()),
      partition(partition),
      kinetic_solver(reactions, partition),
      smart_kinetic_solver(reactions, partition)
    {
    }

    auto setOptions(const KineticPathOptions& _options) -> void
    {
        // Initialise the options of the kinetic path
        options = _options;

        // Set options of kinetic solvers
        kinetic_solver.setOptions(options.kinetics);

        // Set options of smart kinetic solvers
        smart_kinetic_solver.setOptions(options.smart_kinetics);

        // Set options of equilibrium solvers
        options.equilibrium = options.equilibrium;
        options.smart_equilibrium = options.smart_equilibrium;

    }

    auto solve(ChemicalState& state, double t0, double t1, const std::string& units) -> void
    {
        t0 = units::convert(t0, units, "s");
        t1 = units::convert(t1, units, "s");

        kinetic_solver.initialize(state, t0);

        double t = t0;

        // Initialize the output of the equilibrium path calculation
        if(output) output.open();

        // Initialize the plots of the equilibrium path calculation
        for(auto& plot : plots) plot.open();

        while(t < t1)
        {
            // Update the output with current state
            if(output) output.update(state, t);

            // Update the plots with current state
            for(auto& plot : plots) plot.update(state, t);

            // Integrate one time step only
            t = kinetic_solver.step(state, t, t1);
        }

        // Update the output with the final state
        if(output) output.update(state, t1);

        // Update the plots with the final state
        for(auto& plot : plots) plot.update(state, t1);
    }

    auto solve(ChemicalState& state, double t0, double dt, int n, const std::string& units) -> void
    {
        t0 = units::convert(t0, units, "s");
        dt = units::convert(dt, units, "s");

        if(!options.use_smart_kinetic_solver)
            kinetic_solver.initialize(state, t0);
        else
            smart_kinetic_solver.initialize(state, t0);

        double t = t0;

        // Step **: Create KineticPathProfiler to track the timing and results of kinetic path
        KineticPathProfiler profiler;

        // Initialize the output of the equilibrium path calculation
        if(output) output.open();

        // Initialize the plots of the equilibrium path calculation
        for(auto& plot : plots) plot.open();

        // Update the output at initial time
        if(output) output.update(state, t);

        // Update the plots at initial time
        for(auto& plot : plots) plot.update(state, t);

        // Loop through the steps
        unsigned int k = 0;

        while(k < n)
        {
            // Integrate one time step only
            if(!options.use_smart_kinetic_solver)
            {

                // Run kinetic calculations
                t = kinetic_solver.solve(state, t, dt);

                // Update the profiler after every call to step method
                profiler.update(kinetic_solver.result());
            }
            else
            {
                // Run smart kinetic calculations
                t = smart_kinetic_solver.solve(state, t, dt);

                // Update the profiler after every call to step method
                profiler.update(smart_kinetic_solver.result());
            }

            // Update the output with current state
            if(output) output.update(state, t);

            // Update the plots with current state
            for(auto& plot : plots) plot.update(state, t);

            // Increase the counter of the steps
            k++;

        }

        // Step **: Collect the analytics related to reactive transport performance
        auto analysis = profiler.analysis();
        auto rt_results = profiler.results();

        if(options.use_smart_equilibrium_solver)
        {
            // Output characteristics of the smart equilibrium solver (e.g., clusters)
            //kinetic_solver.outputSmartSolverInfo();

            // Output to console time and statistics characterising kinetic solver
            //std::cout << profiler;
        }
        if(options.use_smart_kinetic_solver)
        {
            // Output characteristics of the smart kinetic solver (e.g., clusters)
            //smart_kinetic_solver.outputSmartMethodInfo();

            // Output to console time and statistics characterising kinetic solver
            std::cout << profiler;

            // Generate json output file with collected profiling data
            JsonOutput(output.basename() + "-analysis.json") << analysis;
        }
    }

    auto addSource(const ChemicalState& state, double volumerate, const std::string& units) -> void
    {
        if(!options.use_smart_kinetic_solver)
            kinetic_solver.addSource(state, volumerate, units);
        else
            smart_kinetic_solver.addSource(state, volumerate, units);
    }

    auto addPhaseSink(const std::string& phase, double volumerate, const std::string& units) -> void
    {
        if(!options.use_smart_kinetic_solver)
            kinetic_solver.addPhaseSink(phase, volumerate, units);
        else
            smart_kinetic_solver.addPhaseSink(phase, volumerate, units);
    }

    auto addFluidSink(double volumerate, const std::string& units) -> void
    {
        if(!options.use_smart_kinetic_solver)
            kinetic_solver.addFluidSink(volumerate, units);
        else
            smart_kinetic_solver.addFluidSink(volumerate, units);
    }

    auto addSolidSink(double volumerate, const std::string& units) -> void
    {
        if(!options.use_smart_kinetic_solver)
            kinetic_solver.addSolidSink(volumerate, units);
        else
            smart_kinetic_solver.addSolidSink(volumerate, units);
    }

};

KineticPath::KineticPath(const ReactionSystem& reactions)
{
    RuntimeError("Cannot proceed with KineticPath().",
        "KineticPath() constructor is deprecated. "
        "Use constructor KineticPath(const ReactionSystem&, const Partition&) instead.")
}

KineticPath::KineticPath(const ReactionSystem& reactions, const Partition& partition)
: pimpl(new Impl(reactions, partition))
{}

KineticPath::~KineticPath()
{}

auto KineticPath::operator=(KineticPath other) -> KineticPath&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto KineticPath::setOptions(const KineticPathOptions& options) -> void
{
    pimpl->setOptions(options);
}

auto KineticPath::setPartition(const Partition& partition) -> void
{
    RuntimeError("Cannot proceed with KineticPath::setPartition.",
        "KineticPath::setPartition is deprecated. "
        "Use constructor KineticPath(const ReactionSystem&, const Partition&) instead.")
}

auto KineticPath::addSource(const ChemicalState& state, double volumerate, const std::string& units) -> void
{
    pimpl->addSource(state, volumerate, units);
}

auto KineticPath::addPhaseSink(const std::string& phase, double volumerate, const std::string& units) -> void
{
    pimpl->addPhaseSink(phase, volumerate, units);
}

auto KineticPath::addFluidSink(double volumerate, const std::string& units) -> void
{
    pimpl->addFluidSink(volumerate, units);
}

auto KineticPath::addSolidSink(double volumerate, const std::string& units) -> void
{
    pimpl->addSolidSink(volumerate, units);
}

auto KineticPath::solve(ChemicalState& state, double t0, double t1, const std::string& units) -> void
{
    pimpl->solve(state, t0, t1, units);
}

auto KineticPath::solve(ChemicalState& state, double t0, double dt, int n, const std::string& units) -> void
{
    pimpl->solve(state, t0, dt, n, units);
}

auto KineticPath::output() -> ChemicalOutput
{
    pimpl->output = ChemicalOutput(pimpl->reactions);
    return pimpl->output;
}

auto KineticPath::plot() -> ChemicalPlot
{
    pimpl->plots.emplace_back(ChemicalPlot(pimpl->reactions));
    return pimpl->plots.back();
}

auto KineticPath::plots(unsigned num) -> std::vector<ChemicalPlot>
{
    for(unsigned i = 0; i < num; ++i) plot();
    return pimpl->plots;
}

auto KineticPath::system() const -> const ChemicalSystem&
{
    return pimpl->reactions.system();
}

auto KineticPath::reactions() const -> const ReactionSystem&
{
    return pimpl->reactions;
}

auto KineticPath::partition() const -> const Partition&
{
    return pimpl->partition;
}

} // namespace Reaktoro
