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

// C++ includes
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/ChemicalFormula.hpp>

namespace Reaktoro {

/// Return a predicate function that checks whether an item has a given symbol.
inline auto withSymbol(const String& symbol)
{
    return [&](auto&& item) { return item.symbol() == symbol; };
}

/// Return a predicate function that checks whether an item has one of the given symbols.
inline auto withSymbols(const Strings& symbols)
{
    return [&](auto&& item) { return contains(symbols, item.symbol()); };
}

/// Return a predicate function that checks whether an item has a given name.
inline auto withName(const String& name)
{
    return [&](auto&& item) { return item.name() == name; };
}

/// Return a predicate function that checks whether an item has one of the given names.
inline auto withNames(const Strings& names)
{
    return [&](auto&& item) { return contains(names, item.name()); };
}

/// Return a predicate function that checks whether an item has a given formula.
inline auto withFormula(const String& formula)
{
    ChemicalFormula cformula(formula);
    return [&](auto&& item) { return item.formula().equivalent(cformula); };
}

/// Return a predicate function that checks whether an item has one of the given formulas.
inline auto withFormulas(const Strings& formulas)
{
    std::vector<ChemicalFormula> cformulas(formulas.begin(), formulas.end());
    return [&](auto&& item) {
        const auto formula = item.formula();
        return containsfn(cformulas, [&](auto&& f) {
            return formula.equivalent(f); });
    };
}

/// Return a predicate function that checks whether an item has a given tag.
inline auto withTag(const String& tag)
{
    return [&](auto&& item) { return contains(item.tags(), tag); };
}

/// Return a predicate function that checks whether an item has given tags.
inline auto withTags(const Strings& tags)
{
    return [&](auto&& item) { return contained(tags, item.tags()); };
}

/// Return a predicate function that checks whether an item does not have a given tag.
inline auto withoutTag(const String& tag)
{
    return [&](auto&& item) { return !contains(item.tags(), tag); };
}

/// Return a predicate function that checks whether an item has given element symbols.
inline auto withElements(const Strings& symbols)
{
    return [&](auto&& item) {
        for(auto&& [element, _] : item.elements())
            if(!contains(symbols, element.symbol()))
                return false;
        return true;
    };
}

/// Return a predicate function that checks whether an item has element symbols in given formulas.
inline auto withElementsOf(const Strings& formulas)
{
    std::vector<ChemicalFormula> cformulas(formulas.begin(), formulas.end());
    Strings symbols;
    for(auto&& formula : cformulas)
        for(auto&& symbol : formula.symbols())
            symbols.push_back(symbol);
    return withElements(unique(symbols));
}

} // namespace Reaktoro
