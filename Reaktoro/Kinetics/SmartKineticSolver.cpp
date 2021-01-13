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

#include "SmartKineticSolver.hpp"

// Reaktoro includes
#include <Reaktoro/Kinetics/SmartKineticOptions.hpp>
#include <Reaktoro/Kinetics/SmartKineticSolverClustering.hpp>
#include <Reaktoro/Kinetics/SmartKineticSolverPriorityQueue.hpp>
#include <Reaktoro/Kinetics/SmartKineticSolverNN.hpp>

namespace Reaktoro {

SmartKineticSolver::SmartKineticSolver(const ReactionSystem& reactions, const Partition& partition) {
    solverptr = std::make_unique<SmartKineticSolverClustering>(reactions, partition);
}

SmartKineticSolver::~SmartKineticSolver()
{
}

auto SmartKineticSolver::setOptions(const SmartKineticOptions& options) -> void
{
    if(!initialized)
    {
        // Switch the tag as smart kinetic solver is initialized
        initialized = true;

        // Initialize smart kinetic solver depending on the chosen options
        switch(options.method)
        {
            case SmartKineticStrategy::NearestNeighbour:
                solverptr = std::make_unique<SmartKineticSolverNN>(solverptr->reactions(), solverptr->partition());
                break;
            case SmartKineticStrategy::PriorityQueue:
                solverptr = std::make_unique<SmartKineticSolverPriorityQueue>(solverptr->reactions(), solverptr->partition());
                break;
            case SmartKineticStrategy::PriorityQueuePrimary:
                solverptr = std::make_unique<SmartKineticSolverPriorityQueue>(solverptr->reactions(), solverptr->partition());
                break;
            case SmartKineticStrategy::Clustering:
                solverptr = std::make_unique<SmartKineticSolverClustering>(solverptr->reactions(), solverptr->partition());
                break;
            case SmartKineticStrategy::ClusteringExtended:
                solverptr = std::make_unique<SmartKineticSolverClusteringExtended>(solverptr->reactions(), solverptr->partition());
                break;
            default:
                solverptr = std::make_unique<SmartKineticSolverClustering>(solverptr->reactions(), solverptr->partition());
                break;
        }
    }
    // Set options of the smart equilibrium options
    solverptr.get()->setOptions(options);

}

auto SmartKineticSolver::initialize(ChemicalState& state, double tstart) -> void
{
    return solverptr->initialize(state, tstart);
}

auto SmartKineticSolver::addSource(ChemicalState state, double volumerate, const std::string& units) -> void
{
    return solverptr->addSource(state, volumerate, units);
}

auto SmartKineticSolver::addPhaseSink(const std::string& phase, double volumerate, const std::string& units) -> void
{
    return solverptr->addPhaseSink(phase, volumerate, units);
}

auto SmartKineticSolver::addFluidSink(double volumerate, const std::string& units) -> void
{
    return solverptr->addFluidSink(volumerate, units);
}

auto SmartKineticSolver::addSolidSink(double volumerate, const std::string& units) -> void
{
    return solverptr->addSolidSink(volumerate, units);
}

auto SmartKineticSolver::solve(ChemicalState& state, double t, double dt) -> double
{
    return solverptr->solve(state, t, dt);
}

auto SmartKineticSolver::solve(ChemicalState& state, double t, double dt, VectorConstRef b) -> void
{
    return solverptr->solve(state, t, dt, b);
}

auto SmartKineticSolver::result() const -> const SmartKineticResult&
{
    return solverptr->result();
}

auto SmartKineticSolver::outputInfo() const -> void
{
    solverptr->outputInfo();
}

} // namespace Reaktoro

