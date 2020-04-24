// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
    ElementalComposition(Map<Element, double> const& elements);

    /// Construct an ElementalComposition object with given element symbols and respective coefficients.
    ElementalComposition(Map<String, double> const& elements);

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

    /// Convert this ElementalComposition object into a Map<Element, double> object.
    operator Map<Element, double>() const;

    /// Convert this ElementalComposition object into a Map<String, double> object.
    operator Map<String, double>() const;

private:
    /// The elements and their coefficients.
    Map<Element, double> m_elements;

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

} // namespace Reaktoro
