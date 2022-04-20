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
    if(species.size() > 0)
    {
        // Get the name of the first species
        auto full_name = species[0].name();
        // Take the substring of the full name till the '_' symbol (PHREEQC convention)
        surface_name = species[0].name().substr(0, full_name.find("_"));
    }

    // Add species parsing the info about the sites
    addSurfaceSpecies(species);
}

auto ComplexationSurface::initializeCharges() -> void
{
    // The number of sorption species on the surface site
    const auto num_species = species_list.size();

    // Charges of the sorption species on the surface complexation site
    z = ArrayXd::Zero(species_list.size());

    // Initialize charges of the sorption species on the surface complexation site
    for(int i = 0; i < num_species; ++i)
        z[i] = species_list[i].charge();
}

// Initialize equivalence numbers (the charge of ionic bond) for the surface complexation species.
auto ComplexationSurface::initializeEquivalentNumbers() -> void
{
    // TODO: figure out the equivalence number for Hfo_sOH, is it 0 or 1.
    // The number of sorption species on the surface
    const auto num_species = species_list.size();

    // Charges of the sorption species on the surface complexation
    ze = ArrayXd::Zero(species_list.size());

    // Initialize charges of the sorption species on the surface complexation site
    for(auto i = 0; i < num_species; ++i)
        ze[i] = exchangerEquivalentsNumber(species_list[i]);

    ze = abs(z);

    std::cout << "ze = " << ze.transpose() << std::endl;
}

/// Return equivalence number of provided species.
auto ComplexationSurface::exchangerEquivalentsNumber(const Species& species) -> real
{
    // Run through the elements of the current species and return the coefficient of the exchanger
    for(auto [element, coeff] : species.elements())
    {
//        std::cout << "element = " << element.symbol() << std::endl;
//        std::cout << "coeff = "<< coeff << std::endl;
        // Loop over the names of the existing sites to find a matching one
        for(auto [key, site] : surface_sites)
        {
//            std::cout << "key = " << key << std::endl;
//            std::cout << "site = " << site.name() << std::endl;
            if(element.symbol() == site.name())
                return coeff;
        }
    }

    // If none of the elements contained in species coincide with the name of the sites
    errorif(true, "Could not get information about the exchanger equivalents number. "
                  "Ensure the surface complexation phase contains correct species")
}

auto ComplexationSurface::clone() const -> ComplexationSurface
{
    ComplexationSurface copy = *this;
    return copy;
}

auto ComplexationSurface::name() const -> String
{
    return surface_name;
}

auto ComplexationSurface::potential() const -> real
{
    return state().psi;
}

auto ComplexationSurface::species(Index idx) const -> const Species&
{
    return species_list[idx];
}

auto ComplexationSurface::species() const -> const SpeciesList&
{
    return species_list;
}

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
auto ComplexationSurface::surfaceChargeDensity(ArrayXrConstRef x, ArrayXrConstRef z) const -> real {
    // Auxiliary constants
    const auto F = faradayConstant;

    return F * (z*x).sum() / specific_surface_area;
}

/// Return the complexation surface charge.
auto ComplexationSurface::surfaceCharge(ArrayXrConstRef x, ArrayXrConstRef z) const -> real
{
    return (z*x).sum();
}

/// Return the specific surface area.
auto ComplexationSurface::specificSurfaceArea() const -> real
{
    return specific_surface_area;
}

/// Return the mass.
auto ComplexationSurface::mass() const -> real
{
    return surface_mass;
}

/// Return the list of surface sites.
auto ComplexationSurface::sites() const -> std::map<std::string, ComplexationSurfaceSite>
{
    return surface_sites;
}

auto ComplexationSurface::state(real T, real P, ArrayXrConstRef x) -> ComplexationSurfaceState
{
    surface_state.T = T;
    surface_state.P = P;
    surface_state.x = x;
    surface_state.As = specificSurfaceArea();
    surface_state.mass = mass();
    surface_state.z = charges();
    surface_state.Z = surfaceCharge(x, surface_state.z);
    surface_state.sigma = surfaceChargeDensity(x, surface_state.z);

    return surface_state;
}

auto ComplexationSurface::state() const -> ComplexationSurfaceState
{
    return surface_state;
}

auto ComplexationSurface::addSurfaceSpecies(const SpeciesList& species) -> ComplexationSurface&
{
    // Initialize surface's species list
    species_list = species;

    std::cout << "Species of ComplexationSurface:" << std::endl;
    for (auto s : species)
        std::cout << s.name() << " ";
    std::cout << std::endl;

    // Initialize the list of surface sites
    auto site = ComplexationSurfaceSite();

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
            site = addSite(surface_name + tag, tag);
            sites_number++;
        }
        // Add the species to the list of sorption species of a site with a tag 'tag'
        surface_sites[tag].addSorptionSpecies(s, index++);
    }

    // Initialize charges and equivalents
    initializeCharges();
    initializeEquivalentNumbers();

    return *this;
}

auto ComplexationSurface::setName(const String& name) -> ComplexationSurface&
{
    surface_name = name;
    return *this;
}

auto ComplexationSurface::setMineral(const String& mineral_name) -> ComplexationSurface&
{
    mineral = Species(mineral_name);
    return *this;
}

// Set the specific surface site surface area (in m2/kg).
auto ComplexationSurface::setSpecificSurfaceArea(double value, String unit) -> ComplexationSurface&
{
    specific_surface_area = units::convert(value, unit, "m2/kg");
    return *this;
}

// Set the mass of the solid (in kg).
auto ComplexationSurface::setMass(double value, String unit) -> ComplexationSurface&
{
    surface_mass = units::convert(value, unit, "kg");
    return *this;
}

auto ComplexationSurface::addSite(const String& site_name, const String& site_tag) -> ComplexationSurfaceSite&
{
    error(site_name.find(site_tag) == std::string::npos,
          "Provided site name " + site_name + " doesn't contain a specified site tag " + site_tag);

    // Check if no site with the tag `site_tag` exist, add the site
    if(surface_sites.find(site_tag) == surface_sites.end())
    {
        // Create a new site with a given site name and tag
        auto site = ComplexationSurfaceSite(site_name, site_tag);
        site.setSurfaceName(surface_name);
        surface_sites[site_tag] = site;
        site_tags.emplace_back(site_tag);
    }
    return surface_sites[site_tag];
}

auto ComplexationSurface::addSite(const ComplexationSurfaceSite& site) -> ComplexationSurfaceSite&
{
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
        // Store the site_tag in the list of the tags
        site_tags.emplace_back(site_tag);
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
        std::cout << "\t ssa          : " << site.specificSurfaceArea() << std::endl;
        std::cout << "\t mass         : " << site.mass() << std::endl;
        std::cout << "\t # of species : " << site.sorptionSpecies().size() << std::endl;
        std::cout << "\t :: index :: Species" << std::endl;

        for(auto i : site.indicesSorptionSpecies())
            std::cout << "\t :: " << i << "     :: " << surface.species()[i].name() << std::endl;

    }
    return out;
}

} // namespace Reaktoro
