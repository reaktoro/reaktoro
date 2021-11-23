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
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Thermodynamics/Surface/SurfaceSite.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProps;
class ChemicalState;
class Phase;
class SurfaceSiteState;
class Surface;

/// The chemical properties of an aqueous phase.
class SurfaceSiteProps
{
public:
    /// Construct an uninitialized SurfaceSiteProps object with given chemical system.
    explicit SurfaceSiteProps(const SurfaceSite& site, const ChemicalSystem& system);

    /// Construct an SurfaceSiteProps object with given chemical state of the system.
    explicit SurfaceSiteProps(const SurfaceSite& site, const ChemicalState& state);

    /// Construct an SurfaceSiteProps object with given chemical properties of the system.
    explicit SurfaceSiteProps(const SurfaceSite& site, const ChemicalProps& props);

    /// Construct a copy of a SurfaceSiteProps object.
    SurfaceSiteProps(const SurfaceSiteProps& other);

    /// Destroy this SurfaceSiteProps object.
    virtual ~SurfaceSiteProps();

    /// Assign a SurfaceSiteProps object to this object.
    auto operator=(SurfaceSiteProps other) -> SurfaceSiteProps&;

    /// Update the surface site properties with given chemical state of the system.
    auto update(const ChemicalState& state) -> void;

    /// Update the surface site properties with given chemical properties of the system.
    auto update(const ChemicalProps& props) -> void;

    /// Return the amount of an element (in moles).
    auto elementAmount(const StringOrIndex& symbol) const -> real;

    /// Return the amounts of the elements (in moles).
    auto elementAmounts() const -> ArrayXr;

    /// Return the amounts of the species on the surface site surface (in moles).
    auto speciesAmounts() const -> ArrayXr;

    /// Return the amounts of an surface site species (in moles).
    auto speciesAmount(const StringOrIndex& name) const -> real;

    /// Return the fraction of the species on the surface site composition (in eq).
    auto speciesFractions() const -> ArrayXr;

    /// Return the fraction of an surface site species.
    auto speciesFraction(const StringOrIndex& name) const -> real;

    /// Return the base-10 logarithms of the species' activity on the surface.
    auto speciesActivitiesLg() const -> ArrayXr;

    /// Return the base-10 logarithm of the species; activity on the surface.
    auto speciesActivityLg(const StringOrIndex& name) const -> real;

    // Return the surface site state.
    auto surfaceSiteState() const -> SurfaceSiteState;

    // Return the complexation surface.
    auto surfaceSite() const -> SurfaceSite;

    // Return the sum of species charges.
    auto charge() const -> real;

    // Return the surface charge density.
    auto sigma(real Z) const -> real;

    /// Return the underlying Phase object for the complexation surface site phase.
    auto phase() const -> const Phase&;

    /// Return the underlying Phase object for the complexation surface site phase.
    auto extra() const -> const Map<String, Any>&;

    /// Output the properties of the exchange phase to a stream.
    auto output(std::ostream& out) const -> void;

    /// Output the properties of the exchange phase to a file.
    auto output(const String& filename) const -> void;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

/// Output an SurfaceSiteProps object to an output stream.
auto operator<<(std::ostream& out, const SurfaceSiteProps& state) -> std::ostream&;

} // namespace Reaktoro
