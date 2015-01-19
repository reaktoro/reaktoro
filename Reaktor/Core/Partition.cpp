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

#include "Partition.hpp"

// C++ includes
#include <set>

// Reaktor includes
#include <Reaktor/Common/SetUtils.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>

namespace Reaktor {

struct Partition::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The universe of indices
    Indices universe;

    /// The indices of the equilibrium species
    Indices indices_equilibrium_species;

    /// The indices of the kinetic species
    Indices indices_kinetic_species;

    /// The indices of the inert species
    Indices indices_inert_species;

    /// The indices of the elements in the equilibrium partition
    Indices indices_equilibrium_elements;

    /// The indices of the elements in the kinetic partition
    Indices indices_kinetic_elements;

    /// The indices of the elements in the inert partition
    Indices indices_inert_elements;

    Impl() {}

    Impl(const ChemicalSystem& system)
    : system(system)
    {
        universe = range(system.species().size());
    }

    auto setEquilibriumSpecies(const Indices& ispecies) -> void
    {
        indices_equilibrium_species = ispecies;
        indices_inert_species       = difference(indices_inert_species, ispecies);
        indices_kinetic_species     = difference(universe, unify(ispecies, indices_inert_species));
        finalise();
    }

    auto setEquilibriumSpecies(const std::vector<std::string>& species) -> void
    {
        setEquilibriumSpecies(indices(species, system.species()));
    }

    auto setKineticSpecies(const Indices& ispecies) -> void
    {
        indices_kinetic_species     = ispecies;
        indices_inert_species       = difference(indices_inert_species, ispecies);
        indices_equilibrium_species = difference(universe, unify(ispecies, indices_inert_species));
        finalise();
    }

    auto setKineticSpecies(const std::vector<std::string>& species) -> void
    {
        setKineticSpecies(indices(species, system.species()));
    }

    auto setInertSpecies(const Indices& ispecies) -> void
    {
        indices_inert_species       = ispecies;
        indices_equilibrium_species = difference(indices_equilibrium_species, ispecies);
        indices_kinetic_species     = difference(indices_kinetic_species, ispecies);
        finalise();
    }

    auto setInertSpecies(const std::vector<std::string>& species) -> void
    {
        setInertSpecies(indices(species, system.species()));
    }

    auto indicesPhasesWithSpecies(const Indices& ispecies) -> Indices
    {
        std::set<Index> iphases;
        for(const Index& idx : ispecies)
            iphases.insert(system.indexPhaseWithSpecies(idx));
        return Indices(iphases.begin(), iphases.end());
    }

    auto indicesElementsInPartition(const Indices& ispecies) -> Indices
    {
        std::set<Index> indices;
        for(const Index& i : ispecies)
        {
            const Indices& ielements = system.indicesElementsInSpecies(i);
            indices.insert(ielements.begin(), ielements.end());
        }
        return Indices(indices.begin(), indices.end());
    }

    auto finalise() -> void
    {
        // Initialise the indices of the equilibrium, kinetic and inert elements
        indices_equilibrium_elements$ = idxElementsInPartition(idx_equilibrium_species$);
        indices_kinetic_elements$     = idxElementsInPartition(idx_kinetic_species$);
        indices_inert_elements$       = idxElementsInPartition(idx_inert_species$);

        // Initialise the names of the equilibrium, kinetic and inert elements
        equilibrium_elements$ = extract(system$.elements(), idx_equilibrium_elements$);
        kinetic_elements$     = extract(system$.elements(), idx_kinetic_elements$);
        inert_elements$       = extract(system$.elements(), idx_inert_elements$);

        // Initialise the indices of the phases that contains equilibrium, kinetic and inert elements
        idx_phases_with_equilibrium_species$ = idxPhasesWithSpecies(idx_equilibrium_species$);
        idx_phases_with_kinetic_species$     = idxPhasesWithSpecies(idx_kinetic_species$);
        idx_phases_with_inert_species$       = idxPhasesWithSpecies(idx_inert_species$);
    }
};

Partition::Partition(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

Partition::Partition(const Partition& other)
: pimpl(new Impl(*other.pimpl))
{}

Partition::~Partition()
{}

auto Partition::operator=(Partition other) -> Partition&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto Partition::setEquilibriumSpecies(const Indices& ispecies) -> void
{
    pimpl->setEquilibriumSpecies(ispecies);
}

auto Partition::setEquilibriumSpecies(const std::vector<std::string>& species) -> void
{
    pimpl->setEquilibriumSpecies(species);
}

auto Partition::setKineticSpecies(const Indices& ispecies) -> void
{
    pimpl->setKineticSpecies(ispecies);
}

auto Partition::setKineticSpecies(const std::vector<std::string>& species) -> void
{
    pimpl->setKineticSpecies(species);
}

auto Partition::setInertSpecies(const Indices& ispecies) -> void
{
    pimpl->setInertSpecies(ispecies);
}

auto Partition::setInertSpecies(const std::vector<std::string>& species) -> void
{
    pimpl->setInertSpecies(species);
}

auto Partition::indicesEquilibriumSpecies() const -> const Indices&
{
    return pimpl->indices_equilibrium_species;
}

auto Partition::indicesKineticSpecies() const -> const Indices&
{
    return pimpl->indices_kinetic_species;
}

auto Partition::indicesInertSpecies() const -> const Indices&
{
    return pimpl->indices_inert_species;
}

auto Partition::indicesEquilibriumElements() const -> const Indices&
{
    return pimpl->indices_equilibrium_elements;
}

auto Partition::indicesKineticElements() const -> const Indices&
{
    return pimpl->indices_kinetic_elements;
}

auto Partition::indicesInertElements() const -> const Indices&
{
    return pimpl->indices_inert_elements;
}

auto Partition::allEquilibrium(const ChemicalSystem& system) -> Partition
{
    Indices iequilibrium = range(system.species().size());
    return Partition(system, iequilibrium, {}, {});
}

auto Partition::allKinetic(const ChemicalSystem& system) -> Partition
{
    Indices ikinetic = range(system.species().size());
    return Partition(system, {}, ikinetic, {});
}

auto Partition::allEquilibriumExcept(const ChemicalSystem& system, const Indices& ikinetic, const Indices& iinert) -> Partition
{
    Indices iequilibrium = range(system.species().size());
    iequilibrium = difference(iequilibrium, unify(ikinetic, iinert));
    return Partition(system, iequilibrium, ikinetic, iinert);
}

auto Partition::allKineticExcept(const ChemicalSystem& system, const Indices& iequilibrium, const Indices& iinert) -> Partition
{
    Indices ikinetic = range(system.species().size());
    ikinetic = difference(ikinetic, unify(iequilibrium, iinert));
    return Partition(system, iequilibrium, ikinetic, iinert);
}

} // namespace Reaktor
