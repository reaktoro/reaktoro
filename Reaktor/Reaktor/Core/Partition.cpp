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
#include <Reaktor/Common/Macros.hpp>
#include <Reaktor/Common/SetUtils.hpp>
#include <Reaktor/Core/Multiphase.hpp>
#include <Reaktor/Core/MultiphaseUtils.hpp>

namespace Reaktor {

struct Partition::Impl
{
	/// The sequence of indices [0..numSpecies]
	const Indices universe;

	/// The map that takes the index of a species and return the indices of its elements
	const std::vector<Indices> map_species_to_elements;

    /// The indices of the equilibrium species
    Indices indices_equilibrium_species;

    /// The indices of the kinetic species
    Indices indices_kinetic_species;

    /// The indices of the inert species
    Indices indices_inert_species;

    /// The indices of the equilibrium elements
    Indices indices_equilibrium_elements;

    /// The indices of the equilibrium elements
    Indices indices_kinetic_elements;

    /// The indices of the equilibrium elements
    Indices indices_inert_elements;

	Impl(const Multiphase& multiphase)
	: universe(range(numSpecies(multiphase))),
	  map_species_to_elements(indexMapSpeciesToElements(multiphase)),
	  indices_equilibrium_species(universe),
	  indices_equilibrium_elements(range(numElements(multiphase)))
	{}

	auto setEquilibriumSpecies(const Indices& indices) -> void
    {
        Assert(not indices.empty() and max(indices) < universe.size(),
            "Cannot set the indices of the equilibrium species "
            "with an out-of-bound index.");

        Assert(unique(indices).size() == indices.size(),
            "Cannot set the indices of the equilibrium species "
            "with a non-unique set of indices.");

        indices_equilibrium_species = indices;
        indices_inert_species       = difference(indices_inert_species, indices_equilibrium_species);
        indices_kinetic_species     = difference(universe, unify(indices_equilibrium_species, indices_inert_species));
        finalise();
    }

    auto setKineticSpecies(const Indices& indices) -> void
    {
        Assert(not indices.empty() and max(indices) < universe.size(),
            "Cannot set the indices of the kinetic species "
            "with an out-of-bound index.");

        Assert(unique(indices).size() == indices.size(),
            "Cannot set the indices of the kinetic species "
            "with a non-unique set of indices.");

        indices_kinetic_species     = indices;
        indices_inert_species       = difference(indices_inert_species, indices_kinetic_species);
        indices_equilibrium_species = difference(universe, unify(indices_kinetic_species, indices_inert_species));
        finalise();
    }

    auto setInertSpecies(const Indices& indices) -> void
    {
        Assert(not indices.empty() and max(indices) < universe.size(),
            "Cannot set the indices of the inert species "
            "with an out-of-bound index.");

        Assert(unique(indices).size() == indices.size(),
            "Cannot set the indices of the kinetic species "
            "with a non-unique set of indices.");

        indices_inert_species       = indices;
        indices_equilibrium_species = difference(indices_equilibrium_species, indices_inert_species);
        indices_kinetic_species     = difference(indices_kinetic_species, indices_inert_species);
        finalise();
    }

    auto finalise() -> void
    {
        indices_equilibrium_elements = indicesElementsInPartition(indices_equilibrium_species);
        indices_kinetic_elements = indicesElementsInPartition(indices_kinetic_species);
        indices_inert_elements = indicesElementsInPartition(indices_inert_species);
    }

    auto indicesElementsInPartition(const Indices& ispecies) -> Indices
	{
    	const auto& map = map_species_to_elements;
		std::set<Index> indices;
		for(const Index& i : ispecies)
			indices.insert(map[i].begin(), map[i].end());
		return Indices(indices.begin(), indices.end());
	}
};

Partition::Partition(const Multiphase& multiphase)
: pimpl(new Impl(multiphase))
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

auto Partition::setEquilibriumSpecies(const Indices& indices) -> Partition&
{
	pimpl->setEquilibriumSpecies(indices);
	return *this;
}

auto Partition::setKineticSpecies(const Indices& indices) -> Partition&
{
	pimpl->setKineticSpecies(indices);
	return *this;
}

auto Partition::setInertSpecies(const Indices& indices) -> Partition&
{
	pimpl->setInertSpecies(indices);
	return *this;
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

} /* namespace Reaktor */
