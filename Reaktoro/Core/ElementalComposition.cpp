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

#include "ElementalComposition.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

ElementalComposition::ElementalComposition()
{}

ElementalComposition::ElementalComposition(std::initializer_list<Pair<Element, double>> const& elements)
: m_elements(elements.begin(), elements.end())
{}

ElementalComposition::ElementalComposition(Map<Element, double> const& elements)
: m_elements(elements)
{}

ElementalComposition::ElementalComposition(Map<String, double> const& elements)
{
    for(auto&& [symbol, coeff] : elements)
        m_elements.emplace(Element(symbol), coeff);
}

auto ElementalComposition::size() const -> Index
{
    return m_elements.size();
}

auto ElementalComposition::symbols() const -> Strings
{
    return vectorize(m_elements, RKT_LAMBDA(pair, pair.first.symbol()));
}

auto ElementalComposition::coefficients() const -> Vec<double>
{
    return vectorize(m_elements, RKT_LAMBDA(pair, pair.second));
}

auto ElementalComposition::coefficient(const String& symbol) const -> double
{
    for(auto&& [element, coeff] : m_elements)
        if(element.symbol() == symbol)
            return coeff;
    return 0.0;
}

auto ElementalComposition::molarMass() const -> double
{
    double molar_mass = 0.0;
    for(auto&& [element, coeff] : m_elements)
        molar_mass += element.molarMass() * coeff;
    return molar_mass;
}

ElementalComposition::operator Map<Element, double>() const
{
    return m_elements;
}

ElementalComposition::operator Map<String, double>() const
{
    Map<String, double> map;
    for(auto&& [element, coeff] : m_elements)
        map[element.symbol()] = coeff;
    return map;
}

} // namespace Reaktoro
