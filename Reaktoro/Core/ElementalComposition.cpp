// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
    for(const auto& [symbol, coeff] : elements)
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
    for(const auto& [element, coeff] : m_elements)
        if(element.symbol() == symbol)
            return coeff;
    return 0.0;
}

auto ElementalComposition::molarMass() const -> double
{
    double molar_mass = 0.0;
    for(const auto& [element, coeff] : m_elements)
        molar_mass += element.molarMass() * coeff;
    return molar_mass;
}

auto ElementalComposition::repr() const -> String
{
    std::stringstream ss;
    auto i = 0;
    for(const auto& [element, coeff] : m_elements)
        ss << (i++ == 0 ? "" : " ") << coeff << ":" << element.symbol();
    return ss.str();
}

ElementalComposition::operator Map<Element, double>() const
{
    return m_elements;
}

ElementalComposition::operator Map<String, double>() const
{
    Map<String, double> map;
    for(const auto& [element, coeff] : m_elements)
        map[element.symbol()] = coeff;
    return map;
}

ElementalComposition::operator String() const
{
    return repr();
}

auto parseElementalFormula(const String& formula) -> Pairs<String, double>
{
    auto words = split(formula);
    Pairs<String, double> pairs;
    pairs.reserve(words.size());
    for(auto const& word : words)
    {
        const auto i = word.find(":");
        const auto coeff = tofloat(word.substr(0, i));
        const auto symbol = word.substr(i + 1);
        const auto j = indexfn(pairs, RKT_LAMBDA(x, x.first == symbol)); // check if symbol is already in pairs
        if(j < pairs.size())
            pairs[j].second += coeff; // if symbol alredy in pairs, increment coeff
        else pairs.push_back({symbol, coeff}); // otherwise, insert symbol and its coefficient
    }
    return pairs;
}

} // namespace Reaktoro
