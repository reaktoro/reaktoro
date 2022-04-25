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

/// The chemical properties of an aqueous phase.
class ComplexationSurfaceProps
{
public:
    /// Construct an uninitialized ComplexationSurfaceProps object with given chemical system.
    explicit ComplexationSurfaceProps(const ComplexationSurface& surface, const ChemicalSystem& system);

    /// Construct an ComplexationSurfaceProps object with given chemical state of the system.
    explicit ComplexationSurfaceProps(const ComplexationSurface& surface, const ChemicalState& state);

    /// Construct an ComplexationSurfaceProps object with given chemical properties of the system.
    explicit ComplexationSurfaceProps(const ComplexationSurface& surface, const ChemicalProps& props);

    /// Construct a copy of a ComplexationSurfaceProps object.
    ComplexationSurfaceProps(const ComplexationSurfaceProps& other);

    /// Destroy this ComplexationSurfaceProps object.
    virtual ~ComplexationSurfaceProps();

    /// Assign a ComplexationSurfaceProps object to this object.
    auto operator=(ComplexationSurfaceProps other) -> ComplexationSurfaceProps&;

    /// Update the complexation surface properties with given chemical state of the system.
    auto update(const ChemicalState& state) -> void;

    /// Update the complexation surface properties with given chemical properties of the system.
    auto update(const ChemicalProps& props) -> void;

    /// Return the amount of an element (in moles).
    auto elementAmount(const StringOrIndex& symbol) const -> real;

    /// Return the amounts of the elements (in moles).
    auto elementAmounts() const -> ArrayXr;

    /// Return the amounts of the species on the complexation surface surface (in moles).
    auto speciesAmounts() const -> ArrayXr;

    /// Return the amounts of an complexation surface species (in moles).
    auto speciesAmount(const StringOrIndex& name) const -> real;

    /// Return the fraction of the species on the complexation surface composition (in eq).
    auto speciesFractions() const -> ArrayXr;

    /// Return the fraction of an complexation surface species.
    auto speciesFraction(const StringOrIndex& name) const -> real;

    /// Return the base-10 logarithms of the species' activity on the complexation surface.
    auto speciesActivitiesLg() const -> ArrayXr;

    /// Return the base-10 logarithm of the species; activity on the complexation surface.
    auto speciesActivityLg(const StringOrIndex& name) const -> real;

    // Return the complexation surface state.
    auto complexationSurfaceState() const -> ComplexationSurfaceState;

    // Return the complexation surface.
    auto complexationSurface() const -> ComplexationSurface;

    // Return the surface charge.
    auto surfaceCharge() const -> real;

    // Return the surface charge density.
    auto surfaceChargeDensity(real Z) const -> real;

    // Return the surface potential.
    auto potential(real I, real Z, real sigma) const -> real ;

    /// Return the underlying Phase object for the complexation surface phase.
    auto phase() const -> const Phase&;

    /// Return the underlying Phase object for the complexation surface phase.
    auto extra() const -> const Map<String, Any>&;

    /// Output the properties of the exchange phase to a stream.
    auto output(std::ostream& out) const -> void;

    /// Output the properties of the exchange phase to a file.
    auto output(const String& filename) const -> void;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

/// Output an ComplexationSurfaceProps object to an output stream.
auto operator<<(std::ostream& out, const ComplexationSurfaceProps& state) -> std::ostream&;

} // namespace Reaktoro
