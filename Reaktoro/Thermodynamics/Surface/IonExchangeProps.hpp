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

/// The chemical properties of an aqueous phase.
class IonExchangeProps
{
public:
    /// Construct an uninitialized IonExchangeProps object with given chemical system.
    explicit IonExchangeProps(const ChemicalSystem& system);

    /// Construct an IonExchangeProps object with given chemical state of the system.
    explicit IonExchangeProps(const ChemicalState& state);

    /// Construct an IonExchangeProps object with given chemical properties of the system.
    explicit IonExchangeProps(const ChemicalProps& props);

    /// Construct a copy of a IonExchangeProps object.
    IonExchangeProps(const IonExchangeProps& other);

    /// Destroy this IonExchangeProps object.
    virtual ~IonExchangeProps();

    /// Assign a IonExchangeProps object to this object.
    auto operator=(IonExchangeProps other) -> IonExchangeProps&;

    /// Update the ion exchange properties with given chemical state of the system.
    auto update(const ChemicalState& state) -> void;

    /// Update the ion exchange properties with given chemical properties of the system.
    auto update(const ChemicalProps& props) -> void;

    /// Return the amount of an element (in moles).
    auto elementAmount(const StringOrIndex& symbol) const -> real;

    /// Return the amounts of the elements (in moles).
    auto elementAmounts() const -> ArrayXr;

    /// Return the amounts of the species on the ion exchange surface (in moles).
    auto speciesAmounts() const -> ArrayXr;

    /// Return the amounts of an ion exchange species (in moles).
    auto speciesAmount(const StringOrIndex& name) const -> real;

    /// Return the equivalents of the species on the ion exchange composition (in eq).
    auto speciesEquivalents() const -> ArrayXr;

    /// Return the equivalent of an ion exchange species (in eq).
    auto speciesEquivalent(const StringOrIndex& name) const -> real;

    /// Return the equivalent fractions of the species on the ion exchange surface (in moles) if the molar fractions are provided.
    auto speciesEquivalentFractions() const -> ArrayXr;

    /// Return the equivalent fraction of an ion exchange species.
    auto speciesEquivalentFraction(const StringOrIndex& name) const -> real;

    /// Return the logarithms of the activity coefficients of the species on the ion exchange surface (in moles) if the molar fractions are provided.
    auto speciesActivityCoefficientsLg() const -> ArrayXr;

    /// Return the base-10 logarithm of the activity coefficients of an ion exchange species.
    auto speciesActivityCoefficientLg(const StringOrIndex& name) const -> real;

    /// Return the underlying Phase object for the ion exchange phase.
    auto phase() const -> const Phase&;

    /// Output the properties of the exchange phase to a stream.
    auto output(std::ostream& out) const -> void;

    /// Output the properties of the exchange phase to a file.
    auto output(const String& filename) const -> void;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

/// Output an IonExchangeProps object to an output stream.
auto operator<<(std::ostream& out, const IonExchangeProps& state) -> std::ostream&;

} // namespace Reaktoro
