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
    Impl() {}

    Impl(const Indices& iequilibrium, const Indices& ikinetic, const Indices& iinert)
    : indices_equilibrium_species(iequilibrium),
      indices_kinetic_species(ikinetic),
      indices_inert_species(iinert) {}

    /// The indices of the equilibrium species
    Indices indices_equilibrium_species;

    /// The indices of the kinetic species
    Indices indices_kinetic_species;

    /// The indices of the inert species
    Indices indices_inert_species;
};

Partition::Partition()
: pimpl(new Impl())
{}

Partition::Partition(const Indices& iequilibrium, const Indices& ikinetic, const Indices& iinert)
: pimpl(new Impl(iequilibrium, ikinetic, iinert))
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

auto Partition::equilibriumSpeciesIndices() const -> const Indices&
{
    return pimpl->indices_equilibrium_species;
}

auto Partition::kineticSpeciesIndices() const -> const Indices&
{
    return pimpl->indices_kinetic_species;
}

auto Partition::inertSpeciesIndices() const -> const Indices&
{
    return pimpl->indices_inert_species;
}

auto Partition::allEquilibrium(const Multiphase& multiphase) -> Partition
{
    Indices iequilibrium = range(numSpecies(multiphase));
    return Partition(iequilibrium, Indices(), Indices());
}

auto Partition::allKinetic(const Multiphase& multiphase) -> Partition
{
    Indices ikinetic = range(numSpecies(multiphase));
    return Partition(Indices(), ikinetic, Indices());
}

auto Partition::allEquilibriumExcept(const Multiphase& multiphase, const Indices& ikinetic, const Indices& iinert) -> Partition
{
    Indices iequilibrium = range(numSpecies(multiphase));
    iequilibrium = difference(iequilibrium, unify(ikinetic, iinert));
    return Partition(iequilibrium, ikinetic, iinert);
}

auto Partition::allKineticExcept(const Multiphase& multiphase, const Indices& iequilibrium, const Indices& iinert) -> Partition
{
    Indices ikinetic = range(numSpecies(multiphase));
    ikinetic = difference(ikinetic, unify(iequilibrium, iinert));
    return Partition(iequilibrium, ikinetic, iinert);
}

} // namespace Reaktor
