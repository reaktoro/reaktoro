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

#include "ElementalComposition.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>

namespace Reaktoro {
namespace detail {

/// Return the molar mass of the elemental formula.
auto computeMolarMass(Pairs<Element, double> const& elements) -> double
{
    double molar_mass = 0.0;
    for(auto const& [element, coeff] : elements)
        molar_mass += element.molarMass() * coeff;
    return molar_mass;
}

} // namespace detail

ElementalComposition::ElementalComposition()
{}

ElementalComposition::ElementalComposition(std::initializer_list<Pair<Element, double>> const& elements)
: m_elements(elements.begin(), elements.end()),
  m_molar_mass(detail::computeMolarMass(m_elements))
{}

ElementalComposition::ElementalComposition(Pairs<Element, double> const& elements)
: m_elements(elements),
  m_molar_mass(detail::computeMolarMass(m_elements))
{}

ElementalComposition::ElementalComposition(Pairs<String, double> const& elements)
{
    for(const auto& [symbol, coeff] : elements)
        m_elements.emplace_back(Element(symbol), coeff);
    m_molar_mass = detail::computeMolarMass(m_elements);
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
    for(const auto& [element, coeff] : m_elements)
        if(element.symbol() == symbol)
            return coeff;
    return 0.0;
}

auto ElementalComposition::molarMass() const -> double
{
    return m_molar_mass;
}

auto ElementalComposition::repr() const -> String
{
    std::stringstream ss;
    auto i = 0;
    for(const auto& [element, coeff] : m_elements)
        ss << (i++ == 0 ? "" : " ") << coeff << ":" << element.symbol();
    return ss.str();
}

ElementalComposition::operator Pairs<Element, double>() const
{
    return m_elements;
}

ElementalComposition::operator Pairs<String, double>() const
{
    Pairs<String, double> pairs;
    for(const auto& [element, coeff] : m_elements)
        pairs.emplace_back(element.symbol(), coeff);
    return pairs;
}

ElementalComposition::operator String() const
{
    return repr();
}

auto operator!=(const ElementalComposition& l, const ElementalComposition& r) -> bool
{
    return Map<Element, double>(l) != Map<Element, double>(r);
}

auto operator==(const ElementalComposition& l, const ElementalComposition& r) -> bool
{
    return !(l != r);
}

} // namespace Reaktoro
