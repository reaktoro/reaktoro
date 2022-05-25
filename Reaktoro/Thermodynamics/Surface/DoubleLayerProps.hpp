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
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/ActivityModel.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProps;
class ChemicalState;
class ChemicalSystem;
class Phase;
class Species;
class SpeciesList;

/// The chemical properties of an aqueous phase.
class DoubleLayerProps
{
public:
    /// Construct an uninitialized DoubleLayerProps object with given chemical system.
    explicit DoubleLayerProps(const ChemicalSystem& system);

    /// Construct an DoubleLayerProps object with given chemical state of the system.
    explicit DoubleLayerProps(const ChemicalState& state);

    /// Construct an DoubleLayerProps object with given chemical properties of the system.
    explicit DoubleLayerProps(const ChemicalProps& props);

    /// Construct a copy of a DoubleLayerProps object.
    DoubleLayerProps(const DoubleLayerProps& other);

    /// Destroy this DoubleLayerProps object.
    virtual ~DoubleLayerProps();

    /// Assign a DoubleLayerProps object to this object.
    auto operator=(DoubleLayerProps other) -> DoubleLayerProps&;

    /// Update the double layer properties with given chemical state of the system.
    auto update(const ChemicalState& state) -> void;

    /// Update the double layer properties with given chemical properties of the system.
    auto update(const ChemicalProps& props) -> void;

    /// Return the temperature of the aqueous phase (in K).
    auto temperature() const -> real;

    /// Return the pressure of the aqueous phase (in Pa).
    auto pressure() const -> real;

    /// Return the molality of an element (in molal).
    auto elementMolality(const StringOrIndex& symbol) const -> real;

    /// Return the molality concentrations of the elements in  (in molal).
    auto elementMolalities() const -> ArrayXr;

    /// Return the molality of an double layer solute species (in molal).
    auto speciesMolality(const StringOrIndex& name) const -> real;

    /// Return the molality concentrations of the species (in molal).
    auto speciesMolalities() const -> ArrayXr;

    /// Return the charges of the species.
    auto speciesCharges() const -> ArrayXr;

    /// Return the effective ionic strength of the double layer phase (in molal). Equivalent to @ref ionicStrengthEffective.
    auto ionicStrength() const -> real;

    /// Return the effective ionic strength of the double layer phase (in molal).
    auto ionicStrengthEffective() const -> real;

    /// Return the stoichiometric ionic strength of the double layer phase (in molal).
    auto ionicStrengthStoichiometric() const -> real;

    /// Return the pH of the double layer phase.
    auto pH() const -> real;

    /// Return the pE of the double layer phase.
    auto pE() const -> real;

    /// Return the reduction potential of the double layer phase (in V).
    auto Eh() const -> real;

    /// Return the underlying Phase object for the aqueous phase.
    auto phase() const -> const Phase&;

    /// Output the properties of the aqueous phase to a stream.
    auto output(std::ostream& out) const -> void;

    /// Output the properties of the aqueous phase to a file.
    auto output(const String& filename) const -> void;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

/// Output an DoubleLayerProps object to an output stream.
auto operator<<(std::ostream& out, const DoubleLayerProps& state) -> std::ostream&;

} // namespace Reaktoro
