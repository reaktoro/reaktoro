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

    /// The number of field points
    Index npoints;

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
    Impl(const ChemicalSystem& system, Index npoints)
    : system(system), npoints(npoints),
      states(npoints, ChemicalState(system)),
      equilibriumsolver(system)
    {
        setPartition(Partition(system));
    }

    /// Construct a custom Impl instance with given reaction system
    Impl(const ReactionSystem& reactions, Index npoints)
    : system(reactions.system()), reactions(reactions), npoints(npoints),
      states(npoints, ChemicalState(system)),
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

    auto porosity(ChemicalField& field) const -> void
    {
        field.val.resize(npoints);
        if(field.ddt.rows()) field.ddt.resize(npoints);
        if(field.ddp.rows()) field.ddp.resize(npoints);

        ChemicalScalar porosity;

        for(Index i = 0; i < npoints; ++i)
        {
            porosity = states[i].porosity();

            field.val[i] = porosity.val;
            if(field.ddt.rows()) field.ddt[i] = porosity.ddt;
            if(field.ddp.rows()) field.ddp[i] = porosity.ddp;
        }
    }

    auto saturations(ChemicalField* fields) const -> void
    {
        const Index nfluids = partition.numFluidPhases();

        for(Index j = 0; j < nfluids; ++j)
        {
            fields[j].val.resize(npoints);
            if(fields[j].ddt.rows()) fields[j].ddt.resize(npoints);
            if(fields[j].ddp.rows()) fields[j].ddp.resize(npoints);
        }

        ChemicalVector saturations;

        for(Index i = 0; i < npoints; ++i)
        {
            saturations = states[i].saturations();

            for(Index j = 0; j < nfluids; ++j)
            {
                fields[j].val[i] = saturations[j].val;
                if(fields[j].ddt.rows()) fields[j].ddt[i] = saturations[j].ddt;
                if(fields[j].ddp.rows()) fields[j].ddp[i] = saturations[j].ddp;
            }
        }
    }

    auto densities(ChemicalField* fields) const -> void
    {
        const Index nfluids = partition.numFluidPhases();

        for(Index j = 0; j < nfluids; ++j)
        {
            fields[j].val.resize(npoints);
            if(fields[j].ddt.rows()) fields[j].ddt.resize(npoints);
            if(fields[j].ddp.rows()) fields[j].ddp.resize(npoints);
        }

        ChemicalVector densities;

        for(Index i = 0; i < npoints; ++i)
        {
            densities = states[i].properties().phaseDensities();

            for(Index j = 0; j < nfluids; ++j)
            {
                fields[j].val[i] = densities[j].val;
                if(fields[j].ddt.rows()) fields[j].ddt[i] = densities[j].ddt;
                if(fields[j].ddp.rows()) fields[j].ddp[i] = densities[j].ddp;
            }
        }
    }

    /// Perform one reactive step for every field point.
    auto equilibrate(const double* T, const double* P, const double* be) -> void
    {
        for(Index i = 0; i < states.size(); ++i)
            equilibriumsolver.solve(states[i], T[i], P[i], be + i*Ee);
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

auto ChemicalSolver::porosity(ChemicalField& field) const -> void
{
    return pimpl->porosity(field);
}

auto ChemicalSolver::saturations(ChemicalField* fields) const -> void
{
    return pimpl->saturations(fields);
}

auto ChemicalSolver::densities(ChemicalField* fields) const -> void
{
    return pimpl->densities(fields);
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
