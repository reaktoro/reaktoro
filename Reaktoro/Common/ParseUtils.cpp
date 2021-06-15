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

#include "ParseUtils.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/StringUtils.hpp>

namespace Reaktoro {

auto parseReaction(const String& reaction) -> Pairs<String, double>
{
    Pairs<String, double> equation;
    auto words = split(reaction);
    for(const String& word : words)
    {
        Vec<String> pair = split(word, ":");
        equation.emplace_back(pair[1], tofloat(pair[0]));
    }
    return equation;
}

auto parseNumberStringPairs(const String& str) -> Pairs<String, double>
{
    auto words = split(str);
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
