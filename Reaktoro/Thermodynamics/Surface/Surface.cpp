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

#include "Surface.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcLegacy.hpp>

namespace Reaktoro {

/// ------------------------------------------------------------------------------------------
/// Implementation of the methods of the SurfaceState struct
/// ------------------------------------------------------------------------------------------

/// Update the surface potential for the given ionic strength of the aqueous phase.
auto SurfaceState::updatePotential(real I) -> void
{
    // Using formula sigma = 0.1174*I^0.5*sinh(F*psi/(2*R*T))
    psi = 2*R*T*asinh(sigma/(0.1174*sqrt(I)))/F;
}

/// Update the surface fractions for given indices.
auto SurfaceState::updateFractions(ArrayXrConstRef x_, const Indices& indices) -> void
{
    x(indices) = x_;
}

/// Update the surface charge and charge density.
auto SurfaceState::updateCharge(ArrayXdConstRef z) -> void
{
    charge = (z*x).sum();
    sigma = F*charge/(As*mass);
}

/// ------------------------------------------------------------------------------------------
/// Implementation of the methods of the Surface class
/// ------------------------------------------------------------------------------------------

Surface::Surface()
{}

Surface::Surface(const String& name)
: surface_name(name)
{}

Surface::Surface(const SpeciesList& species)
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

auto Surface::clone() const -> Surface
{
    Surface copy = *this;
    return copy;
}

// Initialize charges for the sorption species.
auto Surface::initializeCharges() -> void
{
    const auto charges = vectorize(species_list, RKT_LAMBDA(x, x.charge()));
    z = ArrayXd::Map(charges.data(), charges.size());
}

// Initialize charges for the sorption site species.
auto Surface::initializeSitesCharges() -> void
{
    for(auto& site : surface_sites) {
        site.second.initializeCharges();
    }
}

// Return the surface name.
auto Surface::name() const -> String
{
    return surface_name;
}

// Return the potential of the surface.
auto Surface::potential() const -> real
{
    return state().psi;
}

// Return the potential of the surface for the given temperature, ionic strength and charge density.
auto Surface::potential(real T, real I, real sigma) const -> real
{
    // Using formula sigma = 0.1174*I^0.5*sinh(F*psi/(2*R*T))
    return 2*R*T*asinh(sigma/(0.1174*sqrt(I)))/F;
}

// Return species of the surface with a given index.
auto Surface::species(Index idx) const -> const Species&
{
    return species_list[idx];
}

// Return species of the surface.
auto Surface::species() const -> const SpeciesList&
{
    return species_list;
}

// Return charges for the surface species.
auto Surface::charges() -> ArrayXd
{
    return z;
}

// Return the mole fractions of the species on surface.
auto Surface::moleFractions() const -> ArrayXr
{
    return surface_state.x;
}

/// Return the surface sigma.
auto Surface::surfaceSigma(real charge) const -> real
{
    return faradayConstant*charge/(ssa*surface_mass);
}

/// Return the surface charge.
auto Surface::surfaceCharge(ArrayXrConstRef x) const -> real
{
    return (z*x).sum();
}

/// Return the specific surface area.
auto Surface::specificSurfaceArea() const -> real
{
    return ssa;
}

/// Return the surface mass.
auto Surface::mass() const -> real
{
    return surface_mass;
}

/// Return the list of surface sites.
auto Surface::sites() const -> std::map<std::string, SurfaceSite>
{
    return surface_sites;
}

/// Return surface state updated for the given temperature, pressure, and zero fractions.
auto Surface::state(real T, real P) -> SurfaceState
{
    surface_state.T = T;
    surface_state.P = P;
    surface_state.x = ArrayXr::Zero(species().size());
    surface_state.As = specificSurfaceArea();
    surface_state.mass = mass();
    surface_state.charge = surfaceCharge(surface_state.x);
    surface_state.sigma = surfaceSigma(surface_state.charge);

    return surface_state;
}

/// Return surface state updated for the given temperature, pressure, and fractions.
auto Surface::state(real T, real P, ArrayXrConstRef x) -> SurfaceState
{
    surface_state.T = T;
    surface_state.P = P;
    surface_state.x = x;
    surface_state.As = specificSurfaceArea();
    surface_state.mass = mass();
    surface_state.charge = surfaceCharge(surface_state.x);
    surface_state.sigma = surfaceSigma(surface_state.charge);

    return surface_state;
}

/// Return surface state.
auto Surface::state() const -> SurfaceState
{
    return surface_state;
}

/// Add list of species forming as the result of sorption on surface.
auto Surface::addSurfaceSpecies(const SpeciesList& species) -> Surface&
{
    // Initialize surface's species list
    species_list = species;

    // Initialize the list of surface sites
    auto site = SurfaceSite();

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

    // Initialize charges of the surface and surface sites
    initializeCharges();
    initializeSitesCharges();

    return *this;
}

/// Set the name of the surface.
auto Surface::setName(const String& name) -> Surface&
{
    surface_name = name;
    return *this;
}

// Set the specific surface site surface area (in m2/kg).
auto Surface::setSpecificSurfaceArea(double value, const String& unit) -> Surface&
{
    ssa = units::convert(value, unit, "m2/kg");
    return *this;
}

// Set the mass of the solid (in kg).
auto Surface::setMass(double value, const String& unit) -> Surface&
{
    surface_mass = units::convert(value, unit, "kg");
    return *this;
}

// Add new site (with a given site name and tag) to the surface.
auto Surface::addSite(const String& site_name, const String& site_tag) -> SurfaceSite&
{
    // Check that provided site name contains the site tag
    error(site_name.find(site_tag) == std::string::npos,
          "Provided site name " + site_name + " doesn't contain a specified site tag " + site_tag);

    // If no site with the tag `site_tag` exist, add the site
    if(surface_sites.find(site_tag) == surface_sites.end())
    {
        // Create a new site with a given site name and tag
        auto site = SurfaceSite(site_name, site_tag);
        site.setSurfaceName(surface_name)
            .setMass(surface_mass)
            .setSpecificSurfaceArea(ssa);

        surface_sites[site_tag] = site;
    }
    return surface_sites[site_tag];
}

// Add the surface site
auto Surface::addSite(const SurfaceSite& site) -> SurfaceSite&
{
    // Fetch the site tag and name from the given site object
    auto site_tag = site.tag();
    auto site_name = site.name();

    // Check if the name of the site is initialized
    error(site_name.empty(), "The name of the site should be initialized.");
    // Check if the name of the site's tag is initialized
    error(site_tag.empty(), "The tag of the site should be initialized.");

    // Check if no site with the tag `site_tag` exist and add either add or extend it with the amount of the site
    auto psite = surface_sites.find(site_tag);
    if(psite == surface_sites.end())
    {
        // Add site into the map with key site_tag
        surface_sites[site_tag] = site;
        // Initialize the name of the surface during the addition of the site
        surface_sites[site_tag].setSurfaceName(surface_name)
                               .setMass(surface_mass)
                               .setSpecificSurfaceArea(ssa);
    }
    else
    {
        // Assign amount to the existing site
        psite->second.setAmount(site.amount());
    }

    return surface_sites[site_tag];
}

/// Output this Surface instance to a stream.
auto Surface::output(std::ostream& out) const -> void
{
    out << *this;
}

/// Output a Surface object to an output stream.
auto operator<<(std::ostream& out, const Surface& surface) -> std::ostream&
{
    std::cout << "Surface: " << surface.name() << std::endl;
    std::cout << "\t specific area, m2/kg : " << surface.specificSurfaceArea() << std::endl;
    std::cout << "\t mass, kg             : " << surface.mass() << std::endl;
    std::cout << "\t # of sites           : " << surface.sites().size() << std::endl;

    for(const auto& [tag, site] : surface.sites())
    {
        std::cout << "Site: " << site.name() << std::endl;
        std::cout << "\t amount       : " << site.amount() << std::endl;
        std::cout << "\t mass         : " << site.mass() << std::endl;
        std::cout << "\t ssa          : " << site.specificSurfaceArea() << std::endl;
        std::cout << "\t # of species : " << site.species().size() << std::endl;
        std::cout << "\t :: index :: species" << std::endl;

        for(auto i : site.speciesIndices())
            std::cout << "\t :: " << i << "     :: " << surface.species()[i].name() << std::endl;
    }
    return out;
}

} // namespace Reaktoro
