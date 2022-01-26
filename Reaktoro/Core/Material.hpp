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
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {

/// A type used to represent a material composed of one or more substances.
class Material
{
public:
    /// Construct a Material object.
    Material(ChemicalSystem const& system);

    /// Set the amount of a substance in the material.
    /// @param substance The chemical formula of the substance.
    /// @param amount The amount of the substance
    /// @param unit The amount unit
    auto setSubstanceAmount(String const& substance, double amount, Chars unit) -> void;

    /// Set the mass of a substance in the material.
    /// @param substance The chemical formula of the substance.
    /// @param mass The mass of the substance
    /// @param unit The mass unit
    auto setSubstanceMass(String const& substance, double mass, Chars unit) -> void;

    /// Add a specified amount or mass of a substance in the material.
    /// @param substance The chemical formula of the substance.
    /// @param value The amount or mass value of the added substance.
    /// @param unit The amount or mass unit (must be convertible to mol or kg).
    /// @warning An error is thrown if `substance` has elements not present in the chemical system.
    auto add(String const& substance, double value, Chars unit) -> void;

    /// Add a specified amount or mass of another material in this material.
    /// @param material The other material composing this material.
    /// @param value The amount or mass value of the added substance.
    /// @param unit The amount or mass unit (must be convertible to mol or kg).
    /// @warning An error is thrown if `material` has elements not present in the chemical system.
    auto add(Material const& material, double value, Chars unit) -> void;

    /// Set a specified amount or mass of a substance in the material.
    /// @param substance The chemical formula of the substance.
    /// @param value The amount or mass value of the added substance.
    /// @param unit The amount or mass unit (must be convertible to mol or kg).
    /// @warning An error is thrown if `substance` has elements not present in the chemical system.
    auto set(String const& substance, double value, Chars unit) -> void;

    /// Set a specified amount or mass of another material in this material.
    /// @param material The other material composing this material.
    /// @param value The amount or mass value of the added substance.
    /// @param unit The amount or mass unit (must be convertible to mol or kg).
    /// @warning An error is thrown if `material` has elements not present in the chemical system.
    auto set(Material const& material, double value, Chars unit) -> void;

    /// Scale the mass of the material to a desired value.
    /// This method will adjust the amounts of the material's
    /// substances to obtain the desired mass value.
    /// @param value The desired mass of the material
    /// @param unit The mass unit
    auto scaleMass(real value, Chars unit) -> void;

    /// Scale the amount of the material to a desired value.
    /// This method will adjust the amounts of the material's
    /// substances to obtain the desired amount value.
    /// @param value The desired amount of the material
    /// @param unit The amount unit
    auto scaleAmount(real value, Chars unit) -> void;

    /// Scale the amount or mass of the material to a desired value.
    /// This method will adjust the amounts of the material's
    /// substances to obtain the desired amount or mass value.
    /// @param value The desired amount of the material
    /// @param unit The amount or mass unit (must be convertible to mol or kg).
    auto scale(real value, Chars unit) -> void;

    /// Return the accumulated amounts of elements and electric charge in the material.
    auto componentAmounts() const -> ArrayXdConstRef;

    /// Return the accumulated amounts of elements in the material.
    auto elementAmounts() const -> ArrayXdConstRef;

    /// Return the accumulated electric charge in the material.
    auto charge() const -> double;

    /// Return the amount of the material as the sum of its substance amounts (in mol).
    auto amount() const -> double;

    /// Return the mass of the material as the sum of its substance masses (in kg).
    auto mass() const -> double;

    // -------------------------------------------------------------------------
    // ARITHMETIC ASSIGNMENT OPERATORS
    // -------------------------------------------------------------------------
    auto operator+=(Material const& other) -> Material&;
    auto operator-=(Material const& other) -> Material&;
    auto operator*=(double scalar) -> Material&;
    auto operator/=(double scalar) -> Material&;

private:
    /// The chemical system associated to this material.
    ChemicalSystem msystem;

    /// The array with the amounts of elements and electric charge.
    ArrayXd b;
};

// -------------------------------------------------------------------------
// ARITHMETIC OPERATORS
// -------------------------------------------------------------------------
auto operator+(const Material& l, const Material& r) -> Material;
auto operator-(const Material& l, const Material& r) -> Material;
auto operator*(double scalar, const Material& r) -> Material;
auto operator*(const Material& l, double scalar) -> Material;
auto operator/(const Material& l, double scalar) -> Material;

} // namespace Reaktoro
