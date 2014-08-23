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
namespace {

template<typename Container>
auto isUnique(Container values) -> bool
{
	std::set<Index> tmp(values.begin(), values.end());
	return tmp.size() == values.size();
}

}

struct Partition::Impl
{
	Impl(const Multiphase& multiphase)
	: universe(range(numSpecies(multiphase))),
	  eindices(universe)
	{}

	/// The sequence of indices [0..numSpecies]
	const Indices universe;

    /// The indices of the equilibrium species
    Indices eindices;

    /// The indices of the kinetic species
    Indices kindices;

    /// The indices of the inert species
    Indices iindices;
};
//
//    auto idxPhasesWithSpecies(const Indices& indices) -> Indices
//    {
//        std::set<Index> idx_phases;
//        for(const Index& idx_species : indices)
//            idx_phases.insert(system$.idxPhaseWithSpecies(idx_species));
//        return Indices(idx_phases.begin(), idx_phases.end());
//    }
//
//    auto idxElementsInPartition(const Indices& idx_species) -> Indices
//    {
//        std::set<Index> indices;
//
//        for(const Index& i : idx_species)
//        {
//            const Indices& idx_elements = system$.idxElementsInSpecies(i);
//            indices.insert(idx_elements.begin(), idx_elements.end());
//        }
//
//        return Indices(indices.begin(), indices.end());
//    }
//
//    auto finalise() -> void
//    {
//        // Initialise the indices of the equilibrium, kinetic and inert elements
//        idx_equilibrium_elements$ = idxElementsInPartition(idx_equilibrium_species$);
//        idx_kinetic_elements$     = idxElementsInPartition(idx_kinetic_species$);
//        idx_inert_elements$       = idxElementsInPartition(idx_inert_species$);
//
//        // Initialise the names of the equilibrium, kinetic and inert elements
//        equilibrium_elements$ = extract(system$.elements(), idx_equilibrium_elements$);
//        kinetic_elements$     = extract(system$.elements(), idx_kinetic_elements$);
//        inert_elements$       = extract(system$.elements(), idx_inert_elements$);
//
//        // Initialise the indices of the phases that contains equilibrium, kinetic and inert elements
//        idx_phases_with_equilibrium_species$ = idxPhasesWithSpecies(idx_equilibrium_species$);
//        idx_phases_with_kinetic_species$     = idxPhasesWithSpecies(idx_kinetic_species$);
//        idx_phases_with_inert_species$       = idxPhasesWithSpecies(idx_inert_species$);
//    }
//
//    auto setEquilibriumSpecies(const std::vector<std::string>& species) -> void
//    {
//        if(species$.empty()) internal::uninitialisedPartitionError();
//
//        internal::checkExistenceOfSpecies(species, species$);
//
//        equilibrium_species$     = species;
//        inert_species$           = difference(inert_species$, equilibrium_species$);
//        kinetic_species$         = difference(species$, unify(equilibrium_species$, inert_species$));
//        idx_equilibrium_species$ = find(equilibrium_species$, species$);
//        idx_kinetic_species$     = find(kinetic_species$, species$);
//        idx_inert_species$       = find(inert_species$, species$);
//
//        finalise();
//    }
//
//    auto setKineticSpecies(const std::vector<std::string>& species) -> void
//    {
//        if(species$.empty()) internal::uninitialisedPartitionError();
//
//        internal::checkExistenceOfSpecies(species, species$);
//
//        kinetic_species$         = species;
//        inert_species$           = difference(inert_species$, kinetic_species$);
//        equilibrium_species$     = difference(species$, unify(kinetic_species$, inert_species$));
//        idx_equilibrium_species$ = find(equilibrium_species$, species$);
//        idx_kinetic_species$     = find(kinetic_species$, species$);
//        idx_inert_species$       = find(inert_species$, species$);
//
//        finalise();
//    }
//
//    auto setInertSpecies(const std::vector<std::string>& species) -> void
//    {
//        if(species$.empty()) internal::uninitialisedPartitionError();
//
//        internal::checkExistenceOfSpecies(species, species$);
//
//        inert_species$           = species;
//        equilibrium_species$     = difference(equilibrium_species$, inert_species$);
//        kinetic_species$         = difference(kinetic_species$, inert_species$);
//        idx_equilibrium_species$ = find(equilibrium_species$, species$);
//        idx_kinetic_species$     = find(kinetic_species$, species$);
//        idx_inert_species$       = find(inert_species$, species$);
//
//        finalise();
//    }
//


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
	const auto& universe = pimpl->universe;
	auto& eindices = pimpl->eindices;
	auto& kindices = pimpl->kindices;
	auto& iindices = pimpl->iindices;

	Assert(not indices.empty() and max(indices) < universe.size(),
		"Cannot set the indices of the equilibrium species "
		"with an out-of-bound index.");

	Assert(unique(indices).size() == indices.size(),
		"Cannot set the indices of the equilibrium species "
		"with a non-unique set of indices.");

	eindices = indices;
	iindices = difference(iindices, eindices);
	kindices = difference(universe, unify(eindices, iindices));

	return *this;
}

auto Partition::setKineticSpecies(const Indices& indices) -> Partition&
{
	const auto& universe = pimpl->universe;
	auto& eindices = pimpl->eindices;
	auto& kindices = pimpl->kindices;
	auto& iindices = pimpl->iindices;

	Assert(not indices.empty() and max(indices) < universe.size(),
		"Cannot set the indices of the kinetic species "
		"with an out-of-bound index.");

	Assert(unique(indices).size() == indices.size(),
		"Cannot set the indices of the kinetic species "
		"with a non-unique set of indices.");

	kindices = indices;
	iindices = difference(iindices, kindices);
	eindices = difference(universe, unify(kindices, iindices));

	return *this;
}

auto Partition::setInertSpecies(const Indices& indices) -> Partition&
{
	const auto& universe = pimpl->universe;
	auto& eindices = pimpl->eindices;
	auto& kindices = pimpl->kindices;
	auto& iindices = pimpl->iindices;

	Assert(not indices.empty() and max(indices) < universe.size(),
		"Cannot set the indices of the inert species "
		"with an out-of-bound index.");

	Assert(unique(indices).size() == indices.size(),
		"Cannot set the indices of the kinetic species "
		"with a non-unique set of indices.");

	iindices = indices;
	eindices = difference(eindices, iindices);
	kindices = difference(kindices, iindices);

	return *this;
}

auto Partition::indicesEquilibriumSpecies() const -> const Indices&
{
	return pimpl->eindices;
}

auto Partition::indicesKineticSpecies() const -> const Indices&
{
	return pimpl->kindices;
}

auto Partition::indicesInertSpecies() const -> const Indices&
{
	return pimpl->iindices;
}

} /* namespace Reaktor */
