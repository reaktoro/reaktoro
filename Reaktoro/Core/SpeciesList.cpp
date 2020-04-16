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

#include "SpeciesList.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>

namespace Reaktoro {
namespace detail {

/// Create ChemicalFormula objects with given formula words.
auto createChemicalFormulas(const Strings& words)
{
    return vectorize(words, lambda(x, ChemicalFormula(x)));
}

} // namespace detail

SpeciesList::SpeciesList()
{}

SpeciesList::SpeciesList(std::initializer_list<Species> species)
: m_species(std::move(species))
{}

SpeciesList::SpeciesList(std::vector<Species> species)
: m_species(std::move(species))
{}

SpeciesList::SpeciesList(StringList formulas)
: m_species(formulas.size())
{
    transform(formulas, m_species, [&](String formula) { return Species(formula); });
}

auto SpeciesList::append(Species substance) -> void
{
    m_species.emplace_back(std::move(substance));
}

auto SpeciesList::data() const -> const std::vector<Species>&
{
    return m_species;
}

auto SpeciesList::size() const -> Index
{
    return data().size();
}

auto SpeciesList::operator[](Index index) const -> const Species&
{
    return data()[index];
}

auto SpeciesList::indexWithName(String name) const -> Index
{
    return indexfn(data(), [&](auto&& s) { return s.name() == name; });
}

auto SpeciesList::indexWithFormula(String substance) const -> Index
{
    ChemicalFormula formula(substance);
    return indexfn(data(), [&](auto&& s) { return s.formula().equivalent(formula); });
}

auto SpeciesList::indexWithSubstance(String substance) const -> Index
{
    return indexfn(data(), [&](auto&& s) { return s.substance() == substance; });
}

auto SpeciesList::withNames(const StringList& names) const -> SpeciesList
{
    return filter(*this, [&](auto&& s) { return contains(names, s.name()); });
}

auto SpeciesList::withFormulas(const StringList& words) const -> SpeciesList
{
    const auto formulas = detail::createChemicalFormulas(words);
    auto pred = [&](auto&& s) {
        const auto formula = s.formula();
        return containsfn(formulas, [&](auto&& f) {
            return formula.equivalent(f); });
    };
    return filter(*this, pred);
}

auto SpeciesList::withSubstances(const StringList& substances) const -> SpeciesList
{
    return filter(*this, [&](auto&& s) { return contains(substances, s.substance()); });
}

auto SpeciesList::withAggregateState(AggregateState state) const -> SpeciesList
{
    return filter(*this, [&](auto&& s) { return s.aggregateState() == state; });
}

auto SpeciesList::withTag(String tag) const -> SpeciesList
{
    return filter(*this, [&](auto&& s) { return contains(s.tags(), tag); });
}


auto SpeciesList::withoutTag(String tag) const -> SpeciesList
{
    return filter(*this, [&](auto&& s) { return !contains(s.tags(), tag); });
}

auto SpeciesList::withTags(const StringList& tags) const -> SpeciesList
{
    return filter(*this, [&](auto&& s) { return contained(tags, s.tags()); });
}

auto SpeciesList::withoutTags(const StringList& tags) const -> SpeciesList
{
    return filter(*this, [&](auto&& s) { return !contained(tags, s.tags()); });
}

auto SpeciesList::withElements(const StringList& symbols) const -> SpeciesList
{
    // Return true if each element symbol in the current species is in symbols
    auto pred = [&](auto&& s) {
        for(auto&& [element, _] : s.elements())
            if(!contains(symbols, element.symbol()))
                return false;
        return true;
    };
    return filter(*this, pred);
}

auto SpeciesList::withElementsOf(const StringList& formulas) const -> SpeciesList
{
    std::vector<String> symbols;
    for(auto&& formula : detail::createChemicalFormulas(formulas))
        for(auto&& symbol : formula.symbols())
            symbols.push_back(symbol);
    return withElements(unique(symbols));
}

} // namespace Reaktoro
