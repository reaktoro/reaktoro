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

#include "SurfaceSite.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Units.hpp>

namespace Reaktoro {

SurfaceSite::SurfaceSite(const String& name)
: site_name(name)
{}

SurfaceSite::SurfaceSite(const String& name, const String& tag)
: site_name(name), site_tag(tag)
{}

// Get the site name.
auto SurfaceSite::name() const -> String
{
    return site_name;
}

// Get the site tag.
auto SurfaceSite::tag() const -> String
{
    return site_tag;
}

// Get the specific surface site surface area (in m2/kg).
auto SurfaceSite::specificSurfaceArea() const -> real
{
    return specific_surface_area;
}

// Get the mass of the surface site (in kg).
auto SurfaceSite::mass() const -> real
{
    return surface_mass;
}

// Get the amount of the surface site (in mol).
auto SurfaceSite::amount() const -> real
{
    return site_amount;
}

// Return sorption surface site species.
auto SurfaceSite::species() const -> const SpeciesList&
{
    return sorption_species;
}

// Return sorption species.
auto SurfaceSite::species(Index idx) const -> const Species&
{
    return sorption_species[idx];
}

// Return indices of the sorption species.
auto SurfaceSite::speciesIndices() const -> Indices
{
    return sorption_species_indices;
}

// Return sorption surface site species' charges.
auto SurfaceSite::charges() const -> ArrayXd
{
    return z;
}

// Initialize charges for the surface site species.
auto SurfaceSite::initializeCharges() -> void
{
    const auto charges = vectorize(sorption_species, RKT_LAMBDA(x, x.charge()));
    z = ArrayXd::Map(charges.data(), charges.size());
}

// Set name of the surface site.
auto SurfaceSite::setName(const String& name) -> SurfaceSite&
{
    site_name = name;

    // If site's name contain `_` symbol,
    // then initialize surface name and site's tag
    auto pos = site_name.find('_');
    if(pos != std::string::npos)
    {
        surface_name = site_name.substr(0, pos);
        site_tag =  site_name.substr(pos, std::string::npos);
    }
    return *this;
}

// Set name of the surface.
auto SurfaceSite::setSurfaceName(const String& name) -> SurfaceSite&
{
    surface_name = name;
    return *this;
}

// Set mass of the surface site.
auto SurfaceSite::setMass(real mass) -> SurfaceSite&
{
    surface_mass = mass;
    return *this;
}

// Set specific surface area of the surface site.
auto SurfaceSite::setSpecificSurfaceArea(real ssa) -> SurfaceSite&
{
    specific_surface_area = ssa;
    return *this;
}

// Set the tag of the surface site.
auto SurfaceSite::setTag(const String& tag) -> SurfaceSite&
{
    // Check if the site name contains this tag
    if(site_name.find(tag) != std::string::npos)
        site_tag = tag;

    return *this;
}

// Set the surface site surface area (in m2).
auto SurfaceSite::setAmount(double value, const String& unit) -> SurfaceSite&
{
    site_amount = units::convert(value, unit, "mol");
    return *this;
}

// Set the density of the solid (in kg).
auto SurfaceSite::setDensity(double value, const String& unit) -> SurfaceSite&
{
    site_density = units::convert(value, unit, "sites/nm2");
    return *this;
}

// Add sorption species on the site.
auto SurfaceSite::addSorptionSpecies(const Species& species, const Index& index) -> void
{
    // If added species' name contains the surface name and the site tag,
    // we add this species and its index to the list of the sorption species and species' indices
    if(species.name().find(surface_name + site_tag) != std::string::npos)
    {
        sorption_species.append(species);
        sorption_species_indices.emplace_back(index);
    }
}

/// Return the surface site sigma.
auto SurfaceSite::siteSigma(real charge) const -> real
{
    return faradayConstant*charge/(specific_surface_area*surface_mass);
}

/// Return the surface site charge.
auto SurfaceSite::siteCharge(ArrayXrConstRef x) -> real
{
     return (z*x).sum();
}

/// Return surface site state updated for the given temperature, pressure, and fractions.
auto SurfaceSite::state(real T, real P, ArrayXrConstRef x) -> SurfaceSiteState
{
    surface_site_state.T = T;
    surface_site_state.P = P;
    surface_site_state.x = x;
    surface_site_state.charge = siteCharge(x);
    surface_site_state.sigma = siteSigma(surface_site_state.charge);

    return surface_site_state;
}

/// Return surface site state.
auto SurfaceSite::state() const -> SurfaceSiteState
{
    return surface_site_state;
}

} // namespace Reaktoro