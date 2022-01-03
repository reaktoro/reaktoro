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
#include <Reaktoro/Core/Element.hpp>

namespace Reaktoro {

/// A type used to describe the elemental composition of chemical species.
class ElementalComposition
{
public:
    /// Construct a default ElementalComposition object.
    ElementalComposition();

    /// Construct an ElementalComposition object with given elements and respective coefficients.
    ElementalComposition(std::initializer_list<Pair<Element, double>> const& elements);

    /// Construct an ElementalComposition object with given elements and respective coefficients.
    ElementalComposition(Pairs<Element, double> const& elements);

    /// Construct an ElementalComposition object with given element symbols and respective coefficients.
    ElementalComposition(Pairs<String, double> const& elements);

    /// Return the number of elements.
    auto size() const -> Index;

    /// Return the symbols of the elements.
    auto symbols() const -> Strings;

    /// Return the coefficients of the elements.
    auto coefficients() const -> Vec<double>;

    /// Return the coefficient of an element symbol in the elemental composition.
    auto coefficient(const String& symbol) const -> double;

    /// Return the molar mass of the elemental composition (in kg/mol).
    auto molarMass() const -> double;

    /// Return a string representation of the elemental composition.
    auto repr() const -> String;

    /// Convert this ElementalComposition object into a Pairs<Element, double> object.
    operator Pairs<Element, double>() const;

    /// Convert this ElementalComposition object into a Pairs<String, double> object.
    operator Pairs<String, double>() const;

    /// Convert this ElementalComposition object into a String representation object.
    operator String() const;

private:
    /// The elements and their coefficients.
    Pairs<Element, double> m_elements;

    /// The molar mass of the elemental formula
    double m_molar_mass = {};

public:
    /// Return begin const iterator of this ElementalComposition instance (for STL compatibility reasons).
    inline auto begin() const { return m_elements.begin(); }

    /// Return begin iterator of this ElementalComposition instance (for STL compatibility reasons).
    inline auto begin() { return m_elements.begin(); }

    /// Return end const iterator of this ElementalComposition instance (for STL compatibility reasons).
    inline auto end() const { return m_elements.end(); }

    /// Return end iterator of this ElementalComposition instance (for STL compatibility reasons).
    inline auto end() { return m_elements.end(); }
};

/// Return true if two ElementalComposition objects are different.
auto operator!=(const ElementalComposition& l, const ElementalComposition& r) -> bool;

/// Return true if two ElementalComposition objects are equal.
auto operator==(const ElementalComposition& l, const ElementalComposition& r) -> bool;

} // namespace Reaktoro
