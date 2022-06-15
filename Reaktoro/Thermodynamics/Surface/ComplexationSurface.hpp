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

// C++ includes
#include <map>

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>
#include <Reaktoro/Thermodynamics/Surface/ComplexationSurfaceSite.hpp>

namespace Reaktoro {

/// A type used to describe the state of the complexation surface.
struct ComplexationSurfaceState
{
    /// Update the surface complexation potential for the given ionic strength of the neighboring phase.
    auto updatePotential(real I) -> void
    {
        // Auxiliary variables
        const auto F = faradayConstant;
        const auto R = universalGasConstant;

        // Using formula sigma = 0.1174*I^0.5*sinh(F*potential/R/T/2) and arcsinh(y) = ln(y+(y^2+1)^1⁄2)
        const auto y = sigma/(0.1174*sqrt(I));
        const auto arcsinhy = std::asinh(y[0]);
        psi = 2*R*T*arcsinhy/F;
    }

    /// The temperature of the solute/gas mixture on the surface (in K).
    real T;

    /// The pressure of the solute/gas mixture on the surface (in Pa).
    real P;

    /// The molar fractions of the complexation surface species.
    ArrayXr x;

    /// The charges of the complexation surface species.
    ArrayXr z;

    /// The complexation surface charge (eq).
    real Z;

    /// The surface charge density (C/m2).
    real sigma;

    // The surface area in (m2/kg).
    real As;

    // The surface mass in (kg).
    real mass;

    // The surface potential (Volt).
    real psi;
};

/// A type used to describe a complexation surface.
/// The ComplexationSurface class is defined as a collection of Species objects, representing,
/// therefore, a composition of complexation phase. Its main purpose is to provide the
/// necessary operations in the calculation of activities of surface complexation species.
class ComplexationSurface
{
    /// Return the complexation surface charge.
    auto surfaceCharge(ArrayXrConstRef x, ArrayXrConstRef z) const -> real;

    // Return the complexation surface charge density.
    auto surfaceChargeDensity(real Z) const -> real;

    // Initialize charges of the surface complexation species.
    auto initializeCharges() -> void;

    // Initialize equivalent numbers of the surface complexation species.
    auto initializeEquivalentNumbers() -> void;

public:

    /// Construct a default ComplexationSurface instance.
    ComplexationSurface();

    /// Construct an ComplexationSurface instance with a given name.
    explicit ComplexationSurface(const String& name);

    /// Construct an ComplexationSurface instance with given species.
    explicit ComplexationSurface(const SpeciesList& species);

    /// Return a deep copy of this ComplexationSurface object.
    auto clone() const -> ComplexationSurface;

    /// Return the name of the surface.
    auto name() const -> String;

    /// Return the value of potential.
    auto potential() const -> real;

    /// Return the species on the surface with given index.
    /// @param idx The index of the species in the complexation surface
    auto species(Index idx) const -> const Species&;

    /// Return the exchange species on the surface.
    auto species() const -> const SpeciesList&;

    /// Return charges of the surface complexation species.
    auto charges() -> ArrayXd;

    /// Return equivalent numbers of surface complexation species.
    auto equivalentsNumbers() -> ArrayXd;

    /// Return the mole fractions of the surface complexation species.
    auto moleFractions() const -> ArrayXr;

    /// Return the specific surface area.
    auto specificSurfaceArea() const -> real;

    /// Return the mass.
    auto mass() const -> real;

    /// Return the list of surface sites.
    auto sites() const -> std::map<std::string, ComplexationSurfaceSite>;

    /// Calculate the state of the aqueous mixture.
    /// @param T The temperature (in K)
    /// @param P The pressure (in Pa)
    /// @param x The fraction of the species in the composition
    auto state(real T, real P, ArrayXrConstRef x) -> ComplexationSurfaceState;

    /// Return the current state of the aqueous mixture.
    auto state() const -> ComplexationSurfaceState;

    /// Add the list species to the surface.
    auto addSurfaceSpecies(const SpeciesList& name) -> ComplexationSurface&;

    /// Set the name of the surface.
    auto setName(const String& surface) -> ComplexationSurface&;

    /// Set the mineral this surface belong to.
    auto setMineral(const String& mineral) -> ComplexationSurface&;

    // Set the specific surface area (in m2/kg).
    auto setSpecificSurfaceArea(double value, const String& unit = "m2/kg") -> ComplexationSurface&;

    // Set the mass of the solid (in kg).
    auto setMass(double value, const String& unit = "kg") -> ComplexationSurface&;

    // Add new site (with a given site name and tag) to the surface.
    auto addSite(const String& site, const String& site_tag) -> ComplexationSurfaceSite&;

    // Add new site to the surface.
    auto addSite(const ComplexationSurfaceSite& site) -> ComplexationSurfaceSite&;

    /// Output this ComplexationSurface instance to a stream.
    auto output(std::ostream& out) const -> void;

private:

    /// The surface complexation state.
    ComplexationSurfaceState surface_state;

    /// Species on the surface complexation surface.
    SpeciesList species_list;

    /// The surface name.
    String surface_name;

    /// The mineral species.
    Species mineral;

    /// The specific area (m2/kg), default value is 600 m2/g = 6e5 m2/kg
    real ssa;

    /// The solid mass (kg)
    real surface_mass;

    /// The charges of the complexation surface species.
    ArrayXr z;

    /// The equivalent numbers of the complexation surface species.
    ArrayXr ze;

    // Map of the surface complexation sites.
    std::map<std::string, ComplexationSurfaceSite> surface_sites;
};

/// Output a ComplexationSurface object to an output stream.
auto operator<<(std::ostream& out, const ComplexationSurface& surface) -> std::ostream&;

} // namespace Reaktoro
