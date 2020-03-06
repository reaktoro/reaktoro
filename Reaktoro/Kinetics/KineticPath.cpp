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

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalOutput.hpp>
#include <Reaktoro/Core/ChemicalPlot.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Kinetics/KineticOptions.hpp>
#include <Reaktoro/Kinetics/KineticSolver.hpp>

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
    KineticSolver kineticsolver;

    /// The options of the kinetic path
    KineticOptions options;

    /// The output instance of the kinetic path calculation
    ChemicalOutput output;

    /// The plots of the kinetic path calculation
    std::vector<ChemicalPlot> plots;

    /// Construct an instance of KineticPath::Impl
    Impl(const ReactionSystem& reactions, const Partition& partition)
    : reactions(reactions),
      system(partition.system()),
      partition(partition),
      kineticsolver(reactions, partition)
    {
    }

    auto setOptions(const KineticOptions& options_) -> void
    {
        // Initialise the options of the kinetic path
        options = options_;

        // Set the options of the kinetic solver
        kineticsolver.setOptions(options);
    }

    auto solve(ChemicalState& state, double t0, double t1, std::string units) -> void
    {
        t0 = units::convert(t0, units, "s");
        t1 = units::convert(t1, units, "s");
        kineticsolver.initialize(state, t0);

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
            t = kineticsolver.step(state, t, t1);
        }

        // Update the output with the final state
        if(output) output.update(state, t1);

        // Update the plots with the final state
        for(auto& plot : plots) plot.update(state, t1);
    }
};

KineticPath::KineticPath(const ReactionSystem& reactions)
{
    RuntimeError("Cannot proceed with KineticPath().",
        "KineticPath() constructor is deprecated. "
        "Use constructor KineticPath(const ReactionSystem&, const Partition&) instead.");
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

auto KineticPath::setOptions(const KineticOptions& options) -> void
{
    pimpl->setOptions(options);
}

auto KineticPath::setPartition(const Partition& partition) -> void
{
    RuntimeError("Cannot proceed with KineticPath::setPartition.",
        "KineticPath::setPartition is deprecated. "
        "Use constructor KineticPath(const ReactionSystem&, const Partition&) instead.");
}

auto KineticPath::addSource(const ChemicalState& state, double volumerate, std::string units) -> void
{
    pimpl->kineticsolver.addSource(state, volumerate, units);
}

auto KineticPath::addPhaseSink(std::string phase, double volumerate, std::string units) -> void
{
    pimpl->kineticsolver.addPhaseSink(phase, volumerate, units);
}

auto KineticPath::addFluidSink(double volumerate, std::string units) -> void
{
    pimpl->kineticsolver.addFluidSink(volumerate, units);
}

auto KineticPath::addSolidSink(double volumerate, std::string units) -> void
{
    pimpl->kineticsolver.addSolidSink(volumerate, units);
}

auto KineticPath::solve(ChemicalState& state, double t0, double t1, std::string units) -> void
{
    pimpl->solve(state, t0, t1, units);
}

auto KineticPath::output() -> ChemicalOutput
{
    pimpl->output = ChemicalOutput(pimpl->reactions);
    return pimpl->output;
}

auto KineticPath::plot() -> ChemicalPlot
{
    pimpl->plots.push_back(ChemicalPlot(pimpl->reactions));
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
