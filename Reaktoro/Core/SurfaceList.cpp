// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

#include "SurfaceList.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {

SurfaceList::SurfaceList()
{}

SurfaceList::SurfaceList(std::initializer_list<Surface> surfaces)
: m_surfaces(std::move(surfaces))
{}

SurfaceList::SurfaceList(Vec<Surface> const& surfaces)
: m_surfaces(surfaces)
{}

auto SurfaceList::append(Surface const& surface) -> void
{
    m_surfaces.push_back(surface);
}

auto SurfaceList::data() const -> Vec<Surface> const&
{
    return m_surfaces;
}

auto SurfaceList::empty() const -> bool
{
    return m_surfaces.empty();
}

auto SurfaceList::size() const -> Index
{
    return m_surfaces.size();
}

auto SurfaceList::operator[](Index i) const -> Surface const&
{
    return m_surfaces[i];
}

auto SurfaceList::operator[](Index i) -> Surface&
{
    return m_surfaces[i];
}

auto SurfaceList::find(String const& symbol) const -> Index
{
    return findWithName(symbol);
}

auto SurfaceList::findWithName(String const& name) const -> Index
{
    return indexfn(m_surfaces, RKT_LAMBDA(e, e.name() == name));
}

auto SurfaceList::index(String const& symbol) const -> Index
{
    return indexWithName(symbol);
}

auto SurfaceList::indexWithName(String const& name) const -> Index
{
    const auto idx = findWithName(name);
    errorif(idx >= size(), "Could not find any Surface object with name ", name, ".");
    return idx;
}

auto SurfaceList::get(String const& symbol) const -> Surface const&
{
    return getWithName(symbol);
}

auto SurfaceList::getWithName(String const& name) const -> Surface const&
{
    return m_surfaces[indexWithName(name)];
}

auto SurfaceList::withNames(StringList const& names) const -> SurfaceList
{
    return vectorize(names, RKT_LAMBDA(name, m_surfaces[indexWithName(name)]));
}

SurfaceList::operator Vec<Surface>&()
{
    return m_surfaces;
}

SurfaceList::operator Vec<Surface> const&() const
{
    return m_surfaces;
}

auto operator+(SurfaceList const& a, SurfaceList const& b) -> SurfaceList
{
    return concatenate(a, b);
}

} // namespace Reaktoro
