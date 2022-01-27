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

#include "SpeciesList.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/ChemicalFormula.hpp>

namespace Reaktoro {

SpeciesList::SpeciesList()
{}

SpeciesList::SpeciesList(std::initializer_list<Species> species)
: m_species(std::move(species))
{}

SpeciesList::SpeciesList(const Vec<Species>& species)
: m_species(species)
{}

SpeciesList::SpeciesList(const StringList& formulas)
: m_species(vectorize(formulas, RKT_LAMBDA(x, Species(x))))
{}

auto SpeciesList::append(const Species& species) -> void
{
    m_species.push_back(species);
}

auto SpeciesList::data() const -> const Vec<Species>&
{
    return m_species;
}

auto SpeciesList::empty() const -> bool
{
    return m_species.empty();
}

auto SpeciesList::size() const -> Index
{
    return m_species.size();
}

auto SpeciesList::elements() const -> ElementList
{
    Map<String, Element> collected;
    for(const auto& s : m_species)
        for(const auto& [element, coeff] : s.elements())
            collected.emplace(element.symbol(), element);
    auto elements = vectorize(collected, RKT_LAMBDA(x, x.second));
    std::sort(elements.begin(), elements.end(),
        [](auto l, auto r) { return l.molarMass() < r.molarMass(); });
    return elements;
}

auto SpeciesList::operator[](Index i) const -> const Species&
{
    return m_species[i];
}

auto SpeciesList::operator[](Index i) -> Species&
{
    return m_species[i];
}

auto SpeciesList::find(const String& name) const -> Index
{
    return findWithName(name);
}

auto SpeciesList::findWithName(const String& name) const -> Index
{
    return indexfn(m_species, RKT_LAMBDA(s, s.name() == name));
}

auto SpeciesList::findWithFormula(const ChemicalFormula& formula) const -> Index
{
    return indexfn(m_species, RKT_LAMBDA(s, formula.equivalent(s.formula())));
}

auto SpeciesList::findWithSubstance(const String& substance) const -> Index
{
    return indexfn(m_species, RKT_LAMBDA(s, s.substance() == substance));
}

auto SpeciesList::index(const String& name) const -> Index
{
    return indexWithName(name);
}

auto SpeciesList::indexWithName(const String& name) const -> Index
{
    const auto idx = findWithName(name);
    error(idx >= size(), "Could not find any Species object with name ", name, ".");
    return idx;
}

auto SpeciesList::indexWithFormula(const ChemicalFormula& formula) const -> Index
{
    const auto idx = findWithFormula(formula);
    error(idx >= size(), "Could not find any Species object with formula ", formula.str(), ".");
    return idx;
}

auto SpeciesList::indexWithSubstance(const String& substance) const -> Index
{
    const auto idx = findWithSubstance(substance);
    error(idx >= size(), "Could not find any Species object with substance ", substance, ".");
    return idx;
}

auto SpeciesList::get(const String& name) const -> const Species&
{
    return getWithName(name);
}

auto SpeciesList::getWithName(const String& name) const -> const Species&
{
    return m_species[indexWithName(name)];
}

auto SpeciesList::getWithFormula(const ChemicalFormula& formula) const -> const Species&
{
    return m_species[indexWithFormula(formula)];
}

auto SpeciesList::getWithSubstance(const String substance) const -> const Species&
{
    return m_species[indexWithSubstance(substance)];
}

auto SpeciesList::withNames(const StringList& names) const -> SpeciesList
{
    return vectorize(names, RKT_LAMBDA(name, m_species[indexWithName(name)]));
}

auto SpeciesList::withFormulas(const StringList& formulas) const -> SpeciesList
{
    return vectorize(formulas, RKT_LAMBDA(formula, m_species[indexWithFormula(formula)]));
}

auto SpeciesList::withSubstances(const StringList& substances) const -> SpeciesList
{
    return vectorize(substances, RKT_LAMBDA(substance, m_species[indexWithSubstance(substance)]));
}

auto SpeciesList::withAggregateState(AggregateState state) const -> SpeciesList
{
    return filter(m_species, RKT_LAMBDA(s, s.aggregateState() == state));
}

auto SpeciesList::withCharge(real value) const -> SpeciesList
{
    return filter(m_species, RKT_LAMBDA(s, s.charge() == value));
}

auto SpeciesList::withTag(String tag) const -> SpeciesList
{
    if(tag.empty())
        return {};
    return filter(m_species, RKT_LAMBDA(s, contains(s.tags(), tag)));
}

auto SpeciesList::withoutTag(String tag) const -> SpeciesList
{
    if(tag.empty())
        return m_species;
    return filter(m_species, RKT_LAMBDA(s, !contains(s.tags(), tag)));
}

auto SpeciesList::withTags(const StringList& tags) const -> SpeciesList
{
    if(tags.empty())
        return {};
    return filter(m_species, RKT_LAMBDA(s, contained(tags, s.tags())));
}

auto SpeciesList::withoutTags(const StringList& tags) const -> SpeciesList
{
    if(tags.empty())
        return m_species;
    return filter(m_species, RKT_LAMBDA(s, !contained(tags, s.tags())));
}

auto SpeciesList::withElements(const StringList& symbols) const -> SpeciesList
{
    return filter(m_species, RKT_LAMBDA(s, contained(s.elements().symbols(), symbols)));
}

auto SpeciesList::withElementsOf(const StringList& formulas) const -> SpeciesList
{
    Strings symbols;
    for(ChemicalFormula formula : formulas)
        symbols = merge(symbols, formula.symbols());
    return withElements(symbols);
}

SpeciesList::operator Vec<Species>&()
{
    return m_species;
}

SpeciesList::operator const Vec<Species>&() const
{
    return m_species;
}

auto operator+(const SpeciesList &a, const SpeciesList &b) -> SpeciesList
{
    return concatenate(a, b);
}

} // namespace Reaktoro
