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
/// @see ComplexationComposition
struct ComplexationSurfaceState
{
    /// Return the surface complexation potential for given ionic strength of the neighboring phase
    auto potential(real I) -> void
    {
        // Auxiliary variables
        const auto F = faradayConstant;
        const auto R = universalGasConstant;

        // Using formula sigma = 0.1174*I^0.5*sinh(F*potential/R/T/2) and arcsinh(y) = ln(y+(y^2+1)^1⁄2)
        const auto y = sigma/(0.1174*sqrt(I));
        psi = 2*R*T/F*log(y + sqrt(1 + y*y));
    }

    /// The temperature of the solute/gas mixture on the surface (in K).
    real T;

    /// The pressure of the solute/gas mixture on the surface (in Pa).
    real P;

    /// The molar fractions of the complexation surface species.
    ArrayXr x;

    /// The charges of the complexation surface species.
    ArrayXr z;

    /// The complexation surface charge (eq) defined as a sum of the species equivalences.
    real Z;

    /// The surface charge density (C/m2) defined as F * Z / As.
    real sigma;

    // The surface area in (m2/kg).
    real As;

    // The surface mass in (kg).
    real mass;

    // The surface potential (Volt).
    real psi = 0.0;
};

/// A type used to describe an complexation surface.
/// The ComplexationSurface class is defined as a collection of Species objects, representing,
/// therefore, a composition of complexation phase. Its main purpose is to provide the
/// necessary operations in the calculation of activities of surface complexation species.
class ComplexationSurface
{
public:

    /// Construct a default ComplexationSurface instance.
    ComplexationSurface();

    /// Construct an ComplexationSurface instance with given name.
    explicit ComplexationSurface(const String& name);

    /// Construct an ComplexationSurface instance with given species.
    explicit ComplexationSurface(const SpeciesList& species);

    /// Return a deep copy of this ComplexationSurface object.
    auto clone() const -> ComplexationSurface;

    /// Return the name of the surface.
    auto name() const -> String;

    /// Return the value of potential.
    auto potential() const -> real;

    /// Return the exchange species on the surface with given index.
    /// @param idx The index of the species in the complexation surface
    auto species(Index idx) const -> const Species&;

    /// Return the exchange species on the surface.
    auto species() const -> const SpeciesList&;

    /// Return charges of the surface complexation species.
    auto charges() -> ArrayXdConstRef;

    /// Return the mole fractions of the surface complexation species.
    auto moleFractions() const -> ArrayXr;

    // Return the complexation surface charge density.
    auto surfaceChargeDensity(ArrayXrConstRef x, ArrayXrConstRef z) const -> real;

    /// Return the complexation surface charge.
    auto surfaceCharge(ArrayXrConstRef x, ArrayXrConstRef z) const -> real;

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

    /// Return the state of the aqueous mixture.
    auto state() const -> ComplexationSurfaceState;

    /// Return the state of the aqueous mixture.
    auto addSurfaceSpecies(const SpeciesList& name) -> ComplexationSurface&;

    /// Set the name of the surface.
    auto setName(const String& surface) -> ComplexationSurface&;

    /// Set the mineral this surface belong to
    auto setMineral(const String& mineral) -> ComplexationSurface&;

    // Add new site (with a given site name and tag) to the surface.
    auto addSite(const String& site, const String& site_tag) -> ComplexationSurfaceSite&;

    // Add new site to the surface.
    auto addSite(const ComplexationSurfaceSite& site) -> ComplexationSurfaceSite&;

    /// Output this ComplexationSurface instance to a stream.
    auto output(std::ostream& out) const -> void;

private:

    /// Current surface complexation state.
    ComplexationSurfaceState surface_state;

    /// All species on the surface complexation surface.
    SpeciesList species_list;

    /// Surface name.
    String surface_name;

    /// Mineral species.
    Species mineral;

    // Map of the surface complexation sites.
    std::map<std::string, ComplexationSurfaceSite> surface_sites;

    // Vector of the surface complexation sites' tags.
    std::vector<std::string> site_tags;

    // The number of sites
    Index sites_number = 0;
};

/// Output a ComplexationSurface object to an output stream.
auto operator<<(std::ostream& out, const ComplexationSurface& surface) -> std::ostream&;

} // namespace Reaktoro
