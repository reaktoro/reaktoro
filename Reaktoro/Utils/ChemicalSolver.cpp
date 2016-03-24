// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

#include "ChemicalSolver.hpp"

// C++ includes
#include <vector>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Kinetics/KineticSolver.hpp>

namespace Reaktoro {

struct ChemicalSolver::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The reaction system instance
    ReactionSystem reactions;

    /// The partitioning of the chemical system
    Partition partition;

    /// The chemical states at each point in the field
    std::vector<ChemicalState> states;

    /// The equilibrium solver
    EquilibriumSolver equilibriumsolver;

    /// The kinetic solver
    KineticSolver kineticsolver;

    /// The porosity field
    ChemicalField phi;

    /// The density field for each fluid phase
    std::vector<ChemicalField> rho;

    /// The saturation field for each fluid phase
    std::vector<ChemicalField> sat;

    /// The sensitivity of the molar amounts of equilibrium species with respect to temperature
    std::vector<Vector> ne_ddt;

    /// The sensitivity of the molar amounts of equilibrium species with respect to pressure
    std::vector<Vector> ne_ddp;

    /// The sensitivity of the molar amounts of equilibrium species with respect to molar amounts of elements
    std::vector<Vector> ne_ddbe;

    /// The number of species and elements in the equilibrium partition
    Index Ne, Ee;

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a custom Impl instance with given chemical system
    Impl(const ChemicalSystem& system, Index size)
    : system(system),
      states(size, ChemicalState(system)),
      equilibriumsolver(system)
    {
        setPartition(Partition(system));
    }

    /// Construct a custom Impl instance with given reaction system
    Impl(const ReactionSystem& reactions, Index size)
    : system(reactions.system()),
      reactions(reactions),
      states(size, ChemicalState(system)),
      equilibriumsolver(system),
      kineticsolver(reactions)
    {
        setPartition(Partition(system));
    }

    /// Set the partition of the chemical system
    auto setPartition(const Partition& partition_) -> void
    {
        partition = partition_;
        equilibriumsolver.setPartition(partition);
        kineticsolver.setPartition(partition);
        Ne = partition.numEquilibriumSpecies();
        Ee = partition.numEquilibriumElements();
    }

    auto porosity() const -> const ChemicalField&
    {
    }

    auto porosityWithDiff() const -> const ChemicalField&
    {
    }

    auto saturation(unsigned ifluidphase) const -> const ChemicalField&
    {
    }

    auto saturationWithDiff(unsigned ifluidphase) const -> const ChemicalField&
    {
    }

    auto density(unsigned ifluidphase) const -> const ChemicalField&
    {
        return rho[ifluidphase];
    }

    auto densityWithDiff(unsigned ifluidphase) const -> const ChemicalField&
    {
    }

    /// Perform one reactive step for every field point.
    auto equilibrate(const double* T, const double* P, const double* be) -> void
    {
        for(Index i = 0; i < states.size(); ++i)
        {
            equilibriumsolver.solve(states[i], T[i], P[i], be + i*Ee);

            const auto& ifp = partition.indicesFluidPhases();
            ChemicalProperties properties = states[i].properties();
            ChemicalVector phase_densities = properties.phaseDensities();
            for(Index j = 0; i < partition.numFluidPhases(); ++i)
            {
                rho[j].val[i] = phase_densities.val[ifp[j]];
                rho[j].ddt[i] = phase_densities.ddt[ifp[j]];
            }
        }
    }

    /// Perform reactive step for every field point.
    auto react(double t, double dt) -> void
    {
        for(Index i = 0; i < states.size(); ++i)
            kineticsolver.solve(states[i], t, dt);
    }
};

ChemicalSolver::ChemicalSolver()
: pimpl(new Impl())
{}

ChemicalSolver::ChemicalSolver(const ChemicalSystem& system, unsigned size)
: pimpl(new Impl(system, size))
{}

ChemicalSolver::ChemicalSolver(const ReactionSystem& reactions, unsigned size)
: pimpl(new Impl(reactions, size))
{}

auto ChemicalSolver::setPartition(const Partition& partition) -> void
{
    pimpl->setPartition(partition);
}

auto ChemicalSolver::setState(const ChemicalState& state) -> void
{
    for(Index i = 0; i < pimpl->states.size(); ++i)
        pimpl->states[i] = state;
}

auto ChemicalSolver::setState(const ChemicalState& state, const Indices& indices) -> void
{
    for(Index i : indices)
        pimpl->states[i] = state;
}

auto ChemicalSolver::porosity() const -> const ChemicalField&
{
    return pimpl->porosity();
}

auto ChemicalSolver::porosityWithDiff() const -> const ChemicalField&
{
    return pimpl->porosityWithDiff();
}

auto ChemicalSolver::saturation(unsigned ifluidphase) const -> const ChemicalField&
{
    return pimpl->saturation(ifluidphase);
}

auto ChemicalSolver::saturationWithDiff(unsigned ifluidphase) const -> const ChemicalField&
{
    return pimpl->saturationWithDiff(ifluidphase);
}

auto ChemicalSolver::density(unsigned ifluidphase) const -> const ChemicalField&
{
    return pimpl->density(ifluidphase);
}

auto ChemicalSolver::densityWithDiff(unsigned ifluidphase) const -> const ChemicalField&
{
    return pimpl->densityWithDiff(ifluidphase);
}

auto ChemicalSolver::equilibrate(const double* T, const double* P, const double* be) -> void
{
    pimpl->equilibrate(T, P, be);
}

auto ChemicalSolver::react(double t, double dt) -> void
{
    pimpl->react(t, dt);
}

}  // namespace Reaktoro
