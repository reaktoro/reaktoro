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
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Kinetics/KineticSolver.hpp>
#include <Reaktoro/Utils/ChemicalField.hpp>

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

    /// The number of species and elements in the equilibrium partition
    Index Ne, Ee, Nk, Nfp;

    /// The number of components
    Index Nc;

    /// The chemical states at each point in the field
    std::vector<ChemicalState> states;

    /// The chemical properties at each point in the field
    std::vector<ChemicalProperties> properties;

    /// The equilibrium solver
    EquilibriumSolver equilibriumsolver;

    /// The kinetic solver
    KineticSolver kineticsolver;

    /// The equilibrium sensitivity at every field point
    std::vector<EquilibriumSensitivity> sensitivities;

    /// The molar amounts of equilibrium species at every field point and their derivatives
    std::vector<ChemicalField> ne;

    /// The porosity at every field point and their derivatives
    ChemicalField porosity;

    /// The saturation of each fluid phase at every field point and their derivatives
    std::vector<ChemicalField> saturations;

    /// The density of each fluid phase at every field point and their derivatives
    std::vector<ChemicalField> densities;

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a custom Impl instance with given chemical system
    Impl(const ChemicalSystem& system, Index npoints)
    : system(system),
      npoints(npoints),
      states(npoints, ChemicalState(system)),
      properties(npoints),
      equilibriumsolver(system)
    {
        setPartition(Partition(system));
    }

    /// Construct a custom Impl instance with given reaction system
    Impl(const ReactionSystem& reactions, Index npoints)
    : system(reactions.system()),
      reactions(reactions),
      npoints(npoints),
      states(npoints, ChemicalState(system)),
      properties(npoints),
      equilibriumsolver(system),
      kineticsolver(reactions)
    {
        setPartition(Partition(system));
    }

    /// Set the partition of the chemical system
    auto setPartition(const Partition& partition_) -> void
    {
        // Set the partition of the chemical solver
        partition = partition_;

        // Initialize the number-type variables
        Ne = partition.numEquilibriumSpecies();
        Nk = partition.numKineticSpecies();
        Ee = partition.numEquilibriumElements();
        Nc = Ee + Nk;
        Nfp = partition.numFluidPhases();

        // Set the partition of the equilibrium and kinetic solvers
        if(Ne) equilibriumsolver.setPartition(partition);
        if(Nk) kineticsolver.setPartition(partition);

        // Initialize the sensitivities member
        sensitivities.resize(npoints);
    }

    /// Update the molar amounts of the equilibrium species and their derivatives at every field point.
    auto updateAmountsEquilibriumSpecies() -> void
    {
        if(ne.size() != Ne)
            ne.resize(Ne, ChemicalField(partition, npoints));

        const Indices& ies = partition.indicesEquilibriumSpecies();

        for(Index k = 0; k < npoints; ++k)
        {
            for(Index i = 0; i < Ne; ++i)
            {
                ne[i].val()[k] = states[k].speciesAmount(ies[i]);
                ne[i].ddT()[k] = sensitivities[k].dnedT[i];
                ne[i].ddP()[k] = sensitivities[k].dnedP[i];
                for(Index j = 0; j < Ee; ++j)
                    ne[i].ddbe()[j][k] = sensitivities[k].dnedbe(i, j);
            }
        }
    }

    /// Update the porosity and their derivatives at every field point.
    auto updatePorosity() -> void
    {
        if(!porosity.size())
            porosity = ChemicalField(partition, npoints);

        ChemicalScalar phi;

        for(Index k = 0; k < npoints; ++k)
        {
            phi = 1.0 - properties[k].solidVolume();
            porosity.set(k, phi, sensitivities[k]);
        }
    }

    /// Get the saturation of each fluid phase at every field point.
    auto updateSaturations() -> void
    {
        if(saturations.size() != Nfp)
            saturations.resize(Nfp, ChemicalField(partition, npoints));

        const Indices& ifp = partition.indicesFluidPhases();

        ChemicalVector fluid_volumes;
        ChemicalVector s;

        for(Index k = 0; k < npoints; ++k)
        {
            fluid_volumes = properties[k].phaseVolumes().rows(ifp);
            s = fluid_volumes/sum(fluid_volumes);
            for(Index j = 0; j < Nfp; ++j)
                saturations[j].set(k, s[j], sensitivities[k]);
        }
    }

    /// Get the density of each fluid phase at every field point.
    auto updateDensities() -> void
    {
        if(densities.size() != Nfp)
            densities.resize(Nfp, ChemicalField(partition, npoints));

        ChemicalVector rho;
        for(Index k = 0; k < npoints; ++k)
        {
            rho = properties[k].phaseDensities();
            for(Index j = 0; j < Nfp; ++j)
                densities[j].set(k, rho[j], sensitivities[k]);
        }
    }

    /// Equilibrate the chemical state at every field point.
    auto equilibrate(const double* T, const double* P, const double* be) -> void
    {
        for(Index i = 0; i < npoints; ++i)
        {
            equilibriumsolver.solve(states[i], T[i], P[i], be + i*Ee);
            properties[i] = states[i].properties();
            sensitivities[i] = equilibriumsolver.sensitivity();
        }
    }

    /// React the chemical state at every field point.
    auto react(double t, double dt) -> void
    {
        for(Index i = 0; i < npoints; ++i)
        {
            kineticsolver.solve(states[i], t, dt);
            properties[i] = states[i].properties();
        }
    }
};

ChemicalSolver::ChemicalSolver()
: pimpl(new Impl())
{}

ChemicalSolver::ChemicalSolver(const ChemicalSystem& system, Index npoints)
: pimpl(new Impl(system, npoints))
{}

ChemicalSolver::ChemicalSolver(const ReactionSystem& reactions, Index npoints)
: pimpl(new Impl(reactions, npoints))
{}

auto ChemicalSolver::setPartition(const Partition& partition) -> void
{
    pimpl->setPartition(partition);
}

auto ChemicalSolver::setState(const ChemicalState& state) -> void
{
    for(Index i = 0; i < pimpl->npoints; ++i)
        pimpl->states[i] = state;
}

auto ChemicalSolver::setStates(const ChemicalState& state, const Indices& indices) -> void
{
    for(Index i : indices)
        pimpl->states[i] = state;
}

auto ChemicalSolver::state(Index i) const -> const ChemicalState&
{
    return pimpl->states[i];
}

auto ChemicalSolver::states() const -> const std::vector<ChemicalState>&
{
    return pimpl->states;
}

auto ChemicalSolver::equilibrate(const double* T, const double* P, const double* be) -> void
{
    pimpl->equilibrate(T, P, be);
}

auto ChemicalSolver::react(double t, double dt) -> void
{
    pimpl->react(t, dt);
}

auto ChemicalSolver::ne() -> const std::vector<ChemicalField>&
{
    pimpl->updateAmountsEquilibriumSpecies();
    return pimpl->ne;
}

auto ChemicalSolver::porosity() -> const ChemicalField&
{
    pimpl->updatePorosity();
    return pimpl->porosity;
}

auto ChemicalSolver::saturations() -> const std::vector<ChemicalField>&
{
    pimpl->updateSaturations();
    return pimpl->saturations;
}

auto ChemicalSolver::densities() -> const std::vector<ChemicalField>&
{
    pimpl->updateDensities();
    return pimpl->densities;
}

}  // namespace Reaktoro
