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

namespace Reaktoro {

// Forward declarations
class ChemicalProps;
class ChemicalState;
class ChemicalSystem;
class Phase;
class ComplexationSurfaceState;
class ComplexationSurface;

/// The chemical properties of an surface and surface sites phase.
class SurfaceProps
{
public:
    /// Construct an uninitialized SurfaceProps object with given chemical system.
    explicit SurfaceProps(const ComplexationSurface& surface, const ChemicalSystem& system);

    /// Construct an SurfaceProps object with given chemical state of the system.
    explicit SurfaceProps(const ComplexationSurface& surface, const ChemicalState& state);

    /// Construct an SurfaceProps object with given chemical properties of the system.
    explicit SurfaceProps(const ComplexationSurface& surface, const ChemicalProps& props);

    /// Construct a copy of a SurfaceProps object.
    SurfaceProps(const SurfaceProps& other);

    /// Destroy this SurfaceProps object.
    virtual ~SurfaceProps();

    /// Assign a SurfaceProps object to this object.
    auto operator=(SurfaceProps other) -> SurfaceProps&;

    /// Update the complexation surface properties with given chemical state of the system.
    auto update(const ChemicalState& state) -> void;

    /// Update the complexation surface properties with given chemical properties of the system.
    auto update(const ChemicalProps& props) -> void;

    /// Return the amount of an element (in moles).
    auto elementAmount(const String& site_tag, const StringOrIndex& symbol) const -> real;

    /// Return the amounts of the elements (in moles).
    auto elementAmounts(const String& site_tag) const -> ArrayXr;

    /// Return the amounts of the species on the complexation surface surface (in moles).
    auto speciesAmounts(const String& site_tag) const -> ArrayXr;

    /// Return the amounts of an complexation surface species (in moles).
    auto speciesAmount(const String& site_tag, const StringOrIndex& name) const -> real;

    /// Return the fraction of the species on the complexation surface composition (in eq).
    auto speciesFractions(const String& site_tag) const -> ArrayXr;

    /// Return the fraction of an complexation surface species.
    auto speciesFraction(const String& site_tag, const StringOrIndex& name) const -> real;

    // Return the complexation surface state.
    auto complexationSurfaceState() const -> ComplexationSurfaceState;

    // Return the complexation surface.
    auto complexationSurface() const -> ComplexationSurface;

    // Return the sum of species charges.
    auto Z() const -> real;

    // Return the surface charge.
    auto charge(real Z) const -> real;

    // Return the surface potential.
    auto potential(real I, real sigma) const -> real ;

    /// Return the underlying Phase object for the complexation surface phase.
    auto phase(const String& site_tag) const -> const Phase&;

    /// Return the underlying Phase object for the complexation surface phase.
    auto extra() const -> const Map<String, Any>&;

    /// Output the properties of the sorption phase to a stream.
    auto output(std::ostream& out) const -> void;

    /// Output the properties of the sorption phase to a file.
    auto output(const String& filename) const -> void;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

/// Output an SurfaceProps object to an output stream.
auto operator<<(std::ostream& out, const SurfaceProps& state) -> std::ostream&;

} // namespace Reaktoro
