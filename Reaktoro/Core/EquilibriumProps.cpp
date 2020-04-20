// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

#include "EquilibriumProps.hpp"

// Reaktoro includes
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {

struct EquilibriumProps::Impl
{
    /// The indices of the equilibrium species partitioned as (primary, secondary).
    ArrayXl ips;

    /// The number of primary species among the equilibrium species.
    Index kp = 0;

    /// The chemical potentials of the species in the equilibrium state (in units of J/mol)
    ArrayXr u;

    /// The chemical potentials of the elements in the equilibrium state (in units of J/mol)
    ArrayXr y;

    /// The stabilities of the species in the equilibrium state (in units of J/mol)
    ArrayXr z;

    /// The indices of equilibrium elements whose amounts should be positive, but given amount was less or equal to zero.
    ArrayXl isue;

    /// The indices of equilibrium species that contain one or more strictly unstable elements.
    ArrayXl isus;

    /// Construct a default EquilibriumProps::Impl instance
    Impl(const ChemicalSystem& system)
    : u(system.species().size()),
      y(system.elements().size()),
      z(system.species().size())
    {}
};

EquilibriumProps::EquilibriumProps(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

EquilibriumProps::EquilibriumProps(const EquilibriumProps& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumProps::~EquilibriumProps()
{}

auto EquilibriumProps::operator=(EquilibriumProps other) -> EquilibriumProps&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumProps::setIndicesEquilibriumSpecies(ArrayXlConstRef ips, Index kp) -> void
{
    pimpl->ips = ips;
    pimpl->kp = kp;
}

auto EquilibriumProps::setIndicesStrictlyUnstableElements(ArrayXlConstRef isue) -> void
{
    pimpl->isue = isue;
}

auto EquilibriumProps::setIndicesStrictlyUnstableSpecies(ArrayXlConstRef isus) -> void
{
    pimpl->isus = isus;
}

auto EquilibriumProps::setSpeciesChemicalPotentials(ArrayXrConstRef u) -> void
{
    assert(u.size() == pimpl->u.size());
    pimpl->u = u;
}

auto EquilibriumProps::setElementChemicalPotentials(ArrayXrConstRef y) -> void
{
    assert(y.size() == pimpl->y.size());
    pimpl->y = y;
}

auto EquilibriumProps::setSpeciesStabilities(ArrayXrConstRef z) -> void
{
    assert(z.size() == pimpl->z.size());
    pimpl->z = z;
}

auto EquilibriumProps::numEquilibriumSpecies() const -> Index
{
    return pimpl->ips.size();
}

auto EquilibriumProps::numPrimarySpecies() const -> Index
{
    return pimpl->kp;
}

auto EquilibriumProps::numSecondarySpecies() const -> Index
{
    return numEquilibriumSpecies() - numPrimarySpecies();
}

auto EquilibriumProps::indicesEquilibriumSpecies() const -> ArrayXlConstRef
{
    return pimpl->ips;
}

auto EquilibriumProps::indicesPrimarySpecies() const -> ArrayXlConstRef
{
    return indicesEquilibriumSpecies().head(numPrimarySpecies());
}

auto EquilibriumProps::indicesSecondarySpecies() const -> ArrayXlConstRef
{
    return indicesEquilibriumSpecies().tail(numSecondarySpecies());
}

auto EquilibriumProps::indicesStrictlyUnstableElements() const -> ArrayXlConstRef
{
    return pimpl->isue;
}

auto EquilibriumProps::indicesStrictlyUnstableSpecies() const -> ArrayXlConstRef
{
    return pimpl->isus;
}

auto EquilibriumProps::speciesChemicalPotentials() const -> ArrayXrConstRef
{
    return pimpl->u;
}

auto EquilibriumProps::elementChemicalPotentials() const -> ArrayXrConstRef
{
    return pimpl->y;
}

auto EquilibriumProps::speciesStabilities() const -> ArrayXrConstRef
{
    return pimpl->z;
}

} // namespace Reaktoro
