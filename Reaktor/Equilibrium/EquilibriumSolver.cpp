// Reaktor is a C++ library for computational reaction modelling.
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

#include "EquilibriumSolver.hpp"

// Reaktor includes
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Partition.hpp>
#include <Reaktor/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktor/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktor/Equilibrium/EquilibriumResult.hpp>

namespace Reaktor {

struct EquilibriumSolver::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The partition of the chemical system in equilibrium and kinetic species
    Partition partition;

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a Impl instance
    Impl(const ChemicalSystem& system)
    : system(system), partition(Partition::allEquilibrium(system))
    {}

    /// Construct a Impl instance
    Impl(const ChemicalSystem& system, const Partition& partition)
    : system(system), partition(partition)
    {}
};

EquilibriumSolver::EquilibriumSolver()
: pimpl(new Impl())
{}

EquilibriumSolver::EquilibriumSolver(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

EquilibriumSolver::EquilibriumSolver(const ChemicalSystem& system, const Partition& partition)
: pimpl(new Impl(system, partition))
{}

auto EquilibriumSolver::solve(const EquilibriumProblem& problem, EquilibriumResult& result) -> void
{
    EquilibriumOptions options;
    solve(problem, result, options);
}

auto EquilibriumSolver::solve(const EquilibriumProblem& problem, EquilibriumResult& result, const EquilibriumOptions& options) -> void
{

}

} // namespace Reaktor
