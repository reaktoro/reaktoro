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
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcLegacy.hpp>
#include <Reaktoro/Common/Constants.hpp>

namespace Reaktoro {

ComplexationSurface::ComplexationSurface()
{}

ComplexationSurface::ComplexationSurface(const String& name)
: surface_name(name)
{}

ComplexationSurface::ComplexationSurface(const SpeciesList& species)
: species_list(species)
{
    // Initialize the surface name from the given species list
    if(species.size() > 0)
    {
        // Get the name of the first species
        auto full_name = species[0].name();
        // Take the substring of the full name till the '_' symbol (PHREEQC convention)
        surface_name = species[0].name().substr(full_name.find("_") + 1);
    }

    // Initialize the list of surface sites
    auto site = ComplexationSurfaceSite();
    Index index = 0;
    for(auto s : species)
    {
        if((s.name().find("_s") != std::string::npos) && (s.charge() == 0))
        {
            // Add a new site indicated with a site_tag "_s" (PHREEQC convention)
            site = addSite(s.name(), "_s");
            sites_number ++;
        }
        else if((s.name().find("_w") != std::string::npos) && (s.charge() == 0))
        {
            // Add a new site indicated with a site_tag "_w" (PHREEQC convention)
            site = addSite(s.name(), "_w");
            sites_number ++;
        }
        else if( ((s.name().find("_w") != std::string::npos) && (s.charge() != 0)) ||
                 ((s.name().find("_s") != std::string::npos) && (s.charge() != 0)))
        {
            // Add a sorbed species to the earlier added site
            site.addSorptionSpecies(s, index);
        }
        else if(s.charge() == 0)
        {
            addSite(s.name(), "");
            sites_number = 1;
        }
        index++;
    }
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
    return state().potential;
}

auto ComplexationSurface::species(Index idx) const -> const Species&
{
    return species_list[idx];
}

auto ComplexationSurface::species() const -> const SpeciesList&
{
    return species_list;
}

/// Initialize the array of species' charges.
auto ComplexationSurface::charges() -> ArrayXdConstRef
{
    // The number of sorption species on the surface site
    const auto num_species = species_list.size();

    // Charges of the sorption species on the surface complexation site
    ArrayXd z = ArrayXd::Zero(species_list.size());

    // Initialize charges of the sorption species on the surface complexation site
    for(int i = 0; i < num_species; ++i)
        z[i] = species_list[i].charge();

    return z;
}

// Return the mole fractions of the species on complexation.
auto ComplexationSurface::moleFractions() const -> ArrayXr
{
    return surface_state.x;
}

/// Return the complexation surface charge density.
auto ComplexationSurface::surfaceChargeDensity(ArrayXrConstRef x, ArrayXrConstRef z) const -> real
{
    real sigma = 0.0;
    const auto F = faradayConstant;

    for (const auto& [tag, site] : surface_sites)
        sigma += F * (z(site.indicesSorptionSpecies())*x(site.indicesSorptionSpecies())).sum() / site.specificSurfaceArea();
    return sigma;
}

/// Return the complexation surface charge.
auto ComplexationSurface::surfaceCharge(ArrayXrConstRef x, ArrayXrConstRef z) const -> real
{
    real Z = 0.0;

    for (const auto& [tag, site] : surface_sites)
        Z += (z(site.indicesSorptionSpecies())*x(site.indicesSorptionSpecies())).sum();
    return Z;
}

/// Return the specific surface area.
auto ComplexationSurface::specificSurfaceArea() const -> real
{
    real ssa = 0;
    for(const auto& [tag, site] : surface_sites)
        ssa += site.specificSurfaceArea();
    return ssa;
}

/// Return the mass.
auto ComplexationSurface::mass() const -> real
{
    real mass = 0;
    for(const auto& [tag, site] : surface_sites)
        mass += site.mass();
    return mass;
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
    surface_state.z = charges();
    surface_state.As = specificSurfaceArea();
    surface_state.Z = surfaceCharge(x, surface_state.z);
    surface_state.sigma = surfaceChargeDensity(x, surface_state.z);
    surface_state.mass = mass();

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

    // Initialize the list of surface sites
    auto site = ComplexationSurfaceSite();

    Index index = 0;
    String aux_tag;

    for(auto s : species)
    {
        // Fetch phreeqc species from the `attachedData` field of the species
        const auto phreeqc_species = std::any_cast<const PhreeqcSpecies*>(s.attachedData());

        if((phreeqc_species->primary != nullptr) && (phreeqc_species->name == phreeqc_species->primary->s->name))
        {
            if( (s.name().find("_s") != std::string::npos) && (surface_sites.find("_s") == surface_sites.end()))
            {
                // Add a new site indicated with a site_tag "_s" (PHREEQC convention)
                site = addSite(surface_name + "_s", "_s");
                sites_number ++;
                aux_tag = "_s";
            }
            else if( (s.name().find("_w") != std::string::npos) && (surface_sites.find("_w") == surface_sites.end()) )
            {
                // Add a new site indicated with a site_tag "_w" (PHREEQC convention)
                site = addSite(surface_name + "_w", "_w");
                sites_number ++;
                aux_tag = "_w";
            }
            else if( (s.name().find("_s") == std::string::npos) || (s.name().find("_w") == std::string::npos))
            {
                site = addSite(s.name(), "");
                sites_number = 1;
                aux_tag = "";
            }
            surface_sites[aux_tag].addSorptionSpecies(s, index);
        }
        else
        {
            // Add a sorbed species to the earlier added site
            if( (s.name().find("_w") != std::string::npos) ||
                (s.name().find("_s") != std::string::npos) )
                surface_sites[aux_tag].addSorptionSpecies(s, index);
        }
        index++;
    }

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

auto ComplexationSurface::addSite(const String& site_name, const String& site_tag) -> ComplexationSurfaceSite&
{
    if(surface_sites.find(site_tag) == surface_sites.end())
    {
        auto site = ComplexationSurfaceSite(site_name, site_tag);
        site.setSurfaceName(surface_name);

        // If previous sites were already added copy their specific surface area and mass
        if(!surface_sites.empty())
        {
            auto ssa = surface_sites[site_tags[0]].specificSurfaceArea();
            auto mass = surface_sites[site_tags[0]].mass();
            site.setSpecificSurfaceArea(ssa, "m2/kg");
            site.setMass(mass, "kg");
        }
        surface_sites[site_tag] = site;
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
    // Hfo
    //	  2.383e-05  Surface charge, eq
    //	  8.611e-03  sigma, C/m2
    //	  5.681e-02  psi, V
    //	 -2.211e+00  -F*psi/RT
    //	  1.096e-01  exp(-F*psif/RT)
    //	  6.000e+01  specific area, m2/kg
    //	  2.670e+02  m2 for   4.450e+00 g

    // Auxiliary variables
    const auto F = faradayConstant;
    const auto R = universalGasConstant;

    std::cout << "Surface: " << surface.name() << std::endl;
    std::cout << "\t charge, eq           : " << surface.state().Z << std::endl;
    std::cout << "\t sigma, C/m2          : " << surface.state().sigma << std::endl;
    std::cout << "\t potential, V         : " << surface.state().potential << std::endl;
    std::cout << "\t -F*psi/RT            : " << - F * surface.state().potential / R / surface.state().T << std::endl;
    std::cout << "\t exp(-F*psi/RT)       : " << exp(- F * surface.state().potential / R / surface.state().T) << std::endl;
    std::cout << "\t specific area, m2/kg : " << surface.specificSurfaceArea() << std::endl;
    std::cout << "\t mass, kg             : " << surface.mass() << std::endl;
    std::cout << "\t # of sites           : " << surface.sites().size() << std::endl;

    // 	  2.500e-05  moles
    //	                                   Mole                     Log
    //	Species               Moles    Fraction    Molality    Molality
    //
    //	Hfo_sOH           1.147e-05       0.459   1.147e-05      -4.940
    //	Hfo_sOHCa+2       1.005e-05       0.402   1.005e-05      -4.998
    //	Hfo_sOH2+         1.754e-06       0.070   1.754e-06      -5.756
    //	Hfo_sO-           1.719e-06       0.069   1.719e-06      -5.765
    //	Hfo_sOCd+         9.469e-10       0.000   9.469e-10      -9.024

    for(const auto& [tag, site] : surface.sites())
    {
        std::cout << "site: " << site.name() << std::endl;
        std::cout << "\t amount       : " << site.amount() << std::endl;
        std::cout << "\t # of species : " << site.sorptionSpecies().size() << std::endl;
        std::cout << "\t :: # :: Species    Mole Fractions" << std::endl;
        for(auto i : site.indicesSorptionSpecies())
        {
//            std::cout << "\t :: " << i << " :: "
//            << surface.species()[i].name() << "    "
//            << surface.moleFractions()[i] << std::endl;
            std::cout << "\t :: " << i << " :: "
                      << surface.species()[i].name() << std::endl;
        }
    }
    return out;
}

} // namespace Reaktoro
