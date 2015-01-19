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
#include <Reaktor/Core/Utils.hpp>

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

    /// The indices of the equilibrium elements
    Indices indices_equilibrium_elements;

    /// The indices of the kinetic elements
    Indices indices_kinetic_elements;

    /// The indices of the inert elements
    Indices indices_inert_elements;

    Impl(const ChemicalSystem& system)
    : system(system)
    {
        universe = range(system.species().size());

        setEquilibriumSpecies(universe);
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
        setEquilibriumSpecies(system.indicesSpecies(species));
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
        setKineticSpecies(system.indicesSpecies(species));
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
        setInertSpecies(system.indicesSpecies(species));
    }

    auto finalise() -> void
    {
        indices_equilibrium_elements = system.indicesElementsInSpecies(indices_equilibrium_species);
        indices_kinetic_elements = system.indicesElementsInSpecies(indices_kinetic_species);
        indices_inert_elements = system.indicesElementsInSpecies(indices_inert_species);
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

auto Partition::indicesElementsInEquilibriumSpecies() const -> const Indices&
{
    return pimpl->indices_equilibrium_elements;
}

auto Partition::indicesElementsInKineticSpecies() const -> const Indices&
{
    return pimpl->indices_kinetic_elements;
}

auto Partition::indicesElementsInInertSpecies() const -> const Indices&
{
    return pimpl->indices_inert_elements;
}

} // namespace Reaktor
