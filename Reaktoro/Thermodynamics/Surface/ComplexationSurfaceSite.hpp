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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

/// A type used to describe the site of a complexation surface.
class ComplexationSurfaceSite
{

public:
    /// Default constructor of ComplexationSurfaceSite class
    ComplexationSurfaceSite() = default;

    /// Copy constructor of ComplexationSurfaceSite class
    ComplexationSurfaceSite(const ComplexationSurfaceSite& other) = default;

    /// Copy assigment operator of ComplexationSurfaceSite class
    ComplexationSurfaceSite& operator=(const ComplexationSurfaceSite& other) = default;

    /// Move constructor of ComplexationSurfaceSite class
    ComplexationSurfaceSite(ComplexationSurfaceSite&& other) = default;

    /// Move assigment operator of ComplexationSurfaceSite class
    ComplexationSurfaceSite& operator=(ComplexationSurfaceSite&& other) = default;

    /// Construct an ComplexationSurfaceSite instance with given name.
    ComplexationSurfaceSite(const String& name);

    /// Construct an ComplexationSurfaceSite instance with given name and tag.
    ComplexationSurfaceSite(const String& name, const String& tag);

    // Get the site name.
    auto name() const -> String;

    // Get the site tag.
    auto tag() const -> String;

    // Get the specific surface site surface area (in m2/kg).
    auto specificSurfaceArea() const -> real;

    // Get the mass of the solid (in kg).
    auto mass() const -> real;

    /// Return the amount of the surface site (in mol).
    auto amount() const -> real;

    // Return the sorption species.
    auto sorptionSpecies() const -> SpeciesList;

    // Return the indices of the sorption species.
    auto sorptionSpeciesIndices() const -> Indices;

    // Set the name of the surface site.
    auto setName(String name) -> ComplexationSurfaceSite&;

    // Set name of the surface.
    auto setSurfaceName(String name) -> ComplexationSurfaceSite&;

    // Set the tag of the surface site.
    auto setTag(String tag) -> ComplexationSurfaceSite&;

    // Set the amount of surface site.
    auto setAmount(double value, String unit = "mol") -> ComplexationSurfaceSite&;

    // Set the density of surface site (in sites/m2).
    auto setDensity(double value, String unit = "sites/m2") -> ComplexationSurfaceSite&;

    // Add sorption species on the site.
    auto addSorptionSpecies(const Species& species, const Index& index) -> void;

private:

    /// The site amount (mol).
    real site_amount;

    /// The site density (sites/nm2, 1e9*sites/m2).
    real site_density;

    /// The name of the surface site.
    String site_name;

    /// The name of the surface.
    String surface_name;

    /// The tag of the surface site.
    String site_tag;

    /// The specific area (m2/kg), default value is 600 m2/g = 6e5 m2/kg
    real specific_surface_area;

    /// The solid mass (kg)
    real surface_mass;

    /// The species that can be sorbed on this site
    SpeciesList sorption_species;

    /// The indices of the species that can be sorbed on this site
    Indices sorption_species_indices;

    // The charges of the species that can be sorbed on this site
    ArrayXd z;

};

}