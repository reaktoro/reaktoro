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

#include "ComplexationSurface.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcLegacy.hpp>

namespace Reaktoro {

ComplexationSurface::ComplexationSurface()
{}

ComplexationSurface::ComplexationSurface(const String& name)
: surface_name(name)
{}

ComplexationSurface::ComplexationSurface(const SpeciesList& species)
{
    // Initialize the surface name from the given species list
    if(!species.empty())
    {
        // Get the name of the first species
        auto full_name = species[0].name();
        // Assign a surface name with a substring of the full name till the '_' symbol (PHREEQC convention)
        surface_name = species[0].name().substr(0, full_name.find("_"));
    }

    // Add species parsing the info about the sites
    addSurfaceSpecies(species);
}

auto ComplexationSurface::clone() const -> ComplexationSurface
{
    ComplexationSurface copy = *this;
    return copy;
}

// Initialize charges for the surface complexation species.
auto ComplexationSurface::initializeCharges() -> void
{
    const auto charges = vectorize(species_list, RKT_LAMBDA(x, x.charge()));
    z = ArrayXd::Map(charges.data(), charges.size());
}

// Initialize equivalence numbers (the charge of ionic bond) for the surface complexation species.
auto ComplexationSurface::initializeEquivalentNumbers() -> void
{
    ze = ArrayXd::Ones(species_list.size());
}

// Return the surface complexation name.
auto ComplexationSurface::name() const -> String
{
    return surface_name;
}

// Return the potential of the surface complexation.
auto ComplexationSurface::potential() const -> real
{
    return state().psi;
}

// Return species of the surface complexation with a given index.
auto ComplexationSurface::species(Index idx) const -> const Species&
{
    return species_list[idx];
}

// Return species of the surface complexation.
auto ComplexationSurface::species() const -> const SpeciesList&
{
    return species_list;
}

// Return charges for the surface complexation species.
auto ComplexationSurface::charges() -> ArrayXd
{
    return z;
}

// Return equivalence numbers for the surface complexation species.
auto ComplexationSurface::equivalentsNumbers() -> ArrayXd
{
    return ze;
}

// Return the mole fractions of the species on complexation.
auto ComplexationSurface::moleFractions() const -> ArrayXr
{
    return surface_state.x;
}

/// Return the complexation surface charge density.
auto ComplexationSurface::surfaceChargeDensity(real Z) const -> real
{
    return faradayConstant*Z/(ssa*surface_mass);
}

/// Return the complexation surface charge.
auto ComplexationSurface::surfaceCharge(ArrayXrConstRef x, ArrayXrConstRef z) const -> real
{
    return (z*x).sum();
}

/// Return the specific surface area.
auto ComplexationSurface::specificSurfaceArea() const -> real
{
    return ssa;
}

/// Return the surface mass.
auto ComplexationSurface::mass() const -> real
{
    return surface_mass;
}

/// Return the list of surface sites.
auto ComplexationSurface::sites() const -> std::map<std::string, ComplexationSurfaceSite>
{
    return surface_sites;
}

/// Return complexation surface state updated for the figen temperature, pressure, and fractions.
auto ComplexationSurface::state(real T, real P, ArrayXrConstRef x) -> ComplexationSurfaceState
{
    surface_state.T = T;
    surface_state.P = P;
    surface_state.x = x;
    surface_state.As = specificSurfaceArea();
    surface_state.mass = mass();
    surface_state.z = charges();
    surface_state.Z = surfaceCharge(x, surface_state.z);
    surface_state.sigma = surfaceChargeDensity(surface_state.Z);

    return surface_state;
}

/// Return complexation surface state.
auto ComplexationSurface::state() const -> ComplexationSurfaceState
{
    return surface_state;
}

/// Add list of species forming as the result of sorption on complexation surface.
auto ComplexationSurface::addSurfaceSpecies(const SpeciesList& species) -> ComplexationSurface&
{
    // Initialize surface's species list
    species_list = species;

    // Initialize the list of surface sites
    auto site = ComplexationSurfaceSite();

    // Auxiliary variables
    Index index = 0;
    String tag;

    for(auto s : species)
    {
        // Fetch phreeqc species from the `attachedData` field of the species
        const auto phreeqc_species = std::any_cast<const PhreeqcSpecies*>(s.attachedData());

        // Get the position of strong or weak sites tag
        auto pos_s = s.name().find("_s");
        auto pos_w = s.name().find("_w");

        // Initialize tag based on the name of the species
        if(pos_s != std::string::npos)
            tag = "_s";
        else if(pos_w != std::string::npos)
            tag = "_w";

        // If the species representing a master species of some site
        if((phreeqc_species->primary != nullptr) &&
           (phreeqc_species->name == phreeqc_species->primary->s->name) &&
           (surface_sites.find(tag) == surface_sites.end())) // if the site with 'tag' tag doesn't exist
        {
            // Add a new site indicated with a site_tag "_s" (PHREEQC convention)
            addSite(surface_name + tag, tag);
        }
        // Add the species to the list of sorption species of a site with a tag 'tag'
        surface_sites[tag].addSorptionSpecies(s, index++);
    }

    // Initialize charges and equivalent numbers
    initializeCharges();
    initializeEquivalentNumbers();

    return *this;
}

/// Set the name of the surface.
auto ComplexationSurface::setName(const String& name) -> ComplexationSurface&
{
    surface_name = name;
    return *this;
}

// Set the mineral associated with the complexation surface.
auto ComplexationSurface::setMineral(const String& mineral_name) -> ComplexationSurface&
{
    mineral = Species(mineral_name);
    return *this;
}

// Set the specific surface site surface area (in m2/kg).
auto ComplexationSurface::setSpecificSurfaceArea(double value, const String& unit) -> ComplexationSurface&
{
    ssa = units::convert(value, unit, "m2/kg");
    return *this;
}

// Set the mass of the solid (in kg).
auto ComplexationSurface::setMass(double value, const String& unit) -> ComplexationSurface&
{
    surface_mass = units::convert(value, unit, "kg");
    return *this;
}

// Add new site (with a given site name and tag) to the surface.
auto ComplexationSurface::addSite(const String& site_name, const String& site_tag) -> ComplexationSurfaceSite&
{
    // Check that provided site name contains the site tag
    error(site_name.find(site_tag) == std::string::npos,
          "Provided site name " + site_name + " doesn't contain a specified site tag " + site_tag);

    // If no site with the tag `site_tag` exist, add the site
    if(surface_sites.find(site_tag) == surface_sites.end())
    {
        // Create a new site with a given site name and tag
        auto site = ComplexationSurfaceSite(site_name, site_tag);
        site.setSurfaceName(surface_name);
        surface_sites[site_tag] = site;
    }
    return surface_sites[site_tag];
}

// Add the complexation surface site
auto ComplexationSurface::addSite(const ComplexationSurfaceSite& site) -> ComplexationSurfaceSite&
{
    // Fetch the site tag and name from the given site object
    auto site_tag = site.tag();
    auto site_name = site.name();

    // Check if the name of the site is initialized
    error(site_name.empty(), "The name of the site should be initialized.");
    // Check if the name of the site's tag is initialized
    error(site_tag.empty(), "The tag of the site should be initialized.");

    // Check if no site with the tag `site_tag` exist, add the site
    if(surface_sites.find(site_tag) == surface_sites.end())
    {
        // Add site into the map with key site_tag
        surface_sites[site_tag] = site;
        // Initialize the name of the surface during the addition of the site
        surface_sites[site_tag].setSurfaceName(surface_name);
    }

    return surface_sites[site_tag];
}

/// Output this ComplexationSurface instance to a stream.
auto ComplexationSurface::output(std::ostream& out) const -> void
{
    out << *this;
}

/// Output a ComplexationSurface object to an output stream.
auto operator<<(std::ostream& out, const ComplexationSurface& surface) -> std::ostream&
{
    std::cout << "SURFACE: " << surface.name() << std::endl;
    std::cout << "\t specific area, m2/kg : " << surface.specificSurfaceArea() << std::endl;
    std::cout << "\t mass, kg             : " << surface.mass() << std::endl;
    std::cout << "\t # of sites           : " << surface.sites().size() << std::endl;

    for(const auto& [tag, site] : surface.sites())
    {
        std::cout << "site: " << site.name() << std::endl;
        std::cout << "\t amount       : " << site.amount() << std::endl;
        std::cout << "\t # of species : " << site.sorptionSpecies().size() << std::endl;
        std::cout << "\t :: index :: Species" << std::endl;

        for(auto i : site.sorptionSpeciesIndices())
            std::cout << "\t :: " << i << "     :: " << surface.species()[i].name() << std::endl;
    }
    return out;
}

} // namespace Reaktoro
