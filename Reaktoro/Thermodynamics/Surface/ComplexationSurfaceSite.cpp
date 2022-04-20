// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

#include "ComplexationSurfaceSite.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Units.hpp>

namespace Reaktoro {

ComplexationSurfaceSite::ComplexationSurfaceSite(const String& name)
: site_name(name)
{}

ComplexationSurfaceSite::ComplexationSurfaceSite(const String& name, const String& tag)
: site_name(name), site_tag(tag)
{}

// Get the site name.
auto ComplexationSurfaceSite::name() const -> String
{
    return site_name;
}

// Get the site tag.
auto ComplexationSurfaceSite::tag() const -> String
{
    return site_tag;
}

// Get the specific surface site surface area (in m2/kg).
auto ComplexationSurfaceSite::specificSurfaceArea() const -> real
{
    return specific_surface_area;
}

// Get the mass of the surface site (in kg).
auto ComplexationSurfaceSite::mass() const -> real
{
    return surface_mass;
}

// Get the amount of the surface site (in mol).
auto ComplexationSurfaceSite::amount() const -> real
{
    return site_amount;
}

// Return sorption species.
auto ComplexationSurfaceSite::sorptionSpecies() const -> SpeciesList
{
    return sorption_species;
}

// Return indices of the sorption species.
auto ComplexationSurfaceSite::indicesSorptionSpecies() const -> Indices
{
    return sorption_species_indices;
}

// Set name of the surface site.
auto ComplexationSurfaceSite::setName(String name) -> ComplexationSurfaceSite&
{
    site_name = name;

    // If site's name contain `_` symbol,
    // then initialize surface name and site's tag
    auto pos = site_name.find("_");
    if(pos != std::string::npos)
    {
        surface_name = site_name.substr(0, pos);
        site_tag =  site_name.substr(pos, std::string::npos);
    }
    return *this;
}

// Set name of the surface.
auto ComplexationSurfaceSite::setSurfaceName(String name) -> ComplexationSurfaceSite&
{
    surface_name = name;
    return *this;
}

// Set the tag of the surface site.
auto ComplexationSurfaceSite::setTag(String tag) -> ComplexationSurfaceSite&
{
    // Check if the site name contains this tag
    if(site_name.find(tag) != std::string::npos)
        site_tag = tag;

    return *this;
}

// Set the surface site surface area (in m2).
auto ComplexationSurfaceSite::setAmount(double value, String unit) -> ComplexationSurfaceSite&
{
    site_amount = units::convert(value, unit, "mol");
    return *this;
}

// Set the density of the solid (in kg).
auto ComplexationSurfaceSite::setDensity(double value, String unit) -> ComplexationSurfaceSite&
{
    site_density = units::convert(value, unit, "sites/nm2");
    return *this;
}

// Add sorption species on the site.
auto ComplexationSurfaceSite::addSorptionSpecies(const Species& species, const Index& index) -> void
{
    // If added species' name contains the surface name and the site tag,
    // we add this species and its index to the list of the sorption species and species' indices
    if(species.name().find(surface_name + site_tag) != std::string::npos)
    {
        sorption_species.append(species);
        sorption_species_indices.emplace_back(index);
    }
}

} // namespace Reaktoro