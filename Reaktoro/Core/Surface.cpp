// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

#include "Surface.hpp"

namespace Reaktoro {

struct Surface::Impl
{
    /// The unique name of the surface.
    String name;

    /// The area model of this surface (in m2).
    SurfaceAreaModel area_model;

    /// Construct a default Surface::Impl object.
    Impl()
    {}
};

Surface::Surface()
: pimpl(new Impl())
{}

Surface::Surface(String const& name)
: Surface()
{
    pimpl->name = name;
}

Surface::Surface(String const& name, SurfaceAreaModel const& model)
: Surface()
{
    pimpl->name = name;
    pimpl->area_model = model;
}

auto Surface::clone() const -> Surface
{
    Surface surface;
    *surface.pimpl = *pimpl;
    return surface;
}

auto Surface::withName(String const& name) const -> Surface
{
    Surface copy = clone();
    copy.pimpl->name = name;
    return copy;
}

auto Surface::withAreaModel(SurfaceAreaModel const& model) const -> Surface
{
    Surface copy = clone();
    copy.pimpl->area_model = model;
    return copy;
}

auto Surface::name() const -> String const&
{
    return pimpl->name;
}

auto Surface::areaModel() const -> SurfaceAreaModel const&
{
    return pimpl->area_model;
}

auto Surface::area(ChemicalProps const& props) const -> real
{
    return pimpl->area_model(props);
}

auto operator<(Surface const& lhs, Surface const& rhs) -> bool
{
    return lhs.name() < rhs.name();
}

auto operator==(Surface const& lhs, Surface const& rhs) -> bool
{
    return lhs.name() == rhs.name();
}

} // namespace Reaktoro
