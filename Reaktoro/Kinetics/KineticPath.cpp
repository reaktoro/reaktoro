// Reaktoro is a C++ library for computational reaction modelling.
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

#include "KineticPath.hpp"

// C++ includes
#include <functional>
using namespace std::placeholders;

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Units.hpp>
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
    KineticSolver solver;

    /// The options of the kinetic path
    KineticOptions options;

    /// The plots of the kinetic path
    std::vector<ChemicalPlot> plots;

    Impl(const ReactionSystem& reactions)
    : reactions(reactions), system(reactions.system()), solver(reactions)
    {
        solver.setPartition(Partition(system));
    }

    auto setOptions(const KineticOptions& options_) -> void
    {
        // Initialise the options of the kinetic path
        options = options_;
    }

    auto setPartition(const Partition& partition) -> void
    {
        solver.setPartition(partition);
    }

    auto setPartition(std::string partition) -> void
    {
        solver.setPartition(partition);
    }

    auto solve(ChemicalState state, double t0, double t1, std::string units) -> void
    {
        t0 = units::convert(t0, units, "s");
        t1 = units::convert(t1, units, "s");
        solver.initialize(state, t0);

        double t = t0;

//        for(unsigned i = 0; i < options.plots.size(); ++i)
//            plots[i].open(options.plots[i]);

        while(t < t1)
        {
            solver.step(state, t, t1);

            // Update the plots
            for(auto& plot : plots)
                plot.update(state, t);
        }
    }
};

KineticPath::KineticPath(const ReactionSystem& reactions)
: pimpl(new Impl(reactions))
{}

KineticPath::KineticPath(const KineticPath& other)
: pimpl(new Impl(*other.pimpl))
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
    pimpl->setPartition(partition);
}

auto KineticPath::setPartition(std::string partition) -> void
{
    pimpl->setPartition(partition);
}

auto KineticPath::solve(const ChemicalState& state, double t0, double t1, std::string units) -> void
{
    pimpl->solve(state, t0, t1, units);
}

} // namespace Reaktoro
