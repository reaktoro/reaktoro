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

#include "Surfaces.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/PhaseList.hpp>
#include <Reaktoro/Core/Surface.hpp>

namespace Reaktoro {

GeneralSurface::GeneralSurface()
{}

GeneralSurface::GeneralSurface(String const& name)
{
    setName(name);
}

auto GeneralSurface::setName(String const& name) -> GeneralSurface&
{
    surface_name = name;
    return *this;
}

auto GeneralSurface::setAreaModel(SurfaceAreaModel const& model) -> GeneralSurface&
{
    area_model = model;
    return *this;
}

auto GeneralSurface::set(SurfaceAreaModel const& model) -> GeneralSurface&
{
    return setAreaModel(model);
}

auto GeneralSurface::name() const -> String const&
{
    return surface_name;
}

auto GeneralSurface::areaModel() const -> SurfaceAreaModel const&
{
    return area_model;
}

auto GeneralSurface::operator()(PhaseList const& phases) const -> Surface
{
    errorif(surface_name.empty(), "Converting a GeneralSurface object to a Surface object requires a non-empty surface name. Use method GeneralSurface::setName to resolve this.");
    errorif(!area_model.initialized(), "Converting a GeneralSurface object to a Surface object requires a non-empty surface area model. Use method GeneralSurface::setAreaModel to resolve this.");

    return Surface()
        .withName(surface_name)
        .withAreaModel(area_model);
}

Surfaces::Surfaces()
{}

auto Surfaces::convert(PhaseList const& phases) const -> Vec<Surface>
{
    Vec<Surface> surfaces;

    for(auto const& fn : surface_generators)
    {
        auto rxns = fn(phases);
        surfaces.insert(surfaces.end(), rxns.begin(), rxns.end());
    }

    return surfaces;
}

} // namespace Reaktoro
