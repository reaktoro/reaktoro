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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/ChemicalFormula.hpp>
#include <Reaktoro/Core/Species.hpp>

namespace Reaktoro {

// Forward declaration of SpeciesListBase
template<typename Data>
class SpeciesListBase;

/// The specialized container to deal with a collection of Species objects.
using SpeciesList = SpeciesListBase<Vec<Species>>;

/// The specialized container to deal with a const reference view of a collection of Species objects.
using SpeciesListConstRef = SpeciesListBase<const Vec<Species>&>;

/// A type used as a collection of species.
template<typename Data>
class SpeciesListBase
{
public:
    /// Construct a default SpeciesListBase object.
    SpeciesListBase()
    {}

    /// Construct an SpeciesListBase object with given species.
    SpeciesListBase(std::initializer_list<Species> species)
     : m_species(std::move(species))
    {}

    /// Construct an SpeciesListBase object with given species.
    SpeciesListBase(const Vec<Species>& species)
     : m_species(species)
    {}

    /// Construct an SpeciesListBase object with given species formulas.
    SpeciesListBase(const StringList& formulas)
     : m_species(vectorize(formulas, RKT_LAMBDA(x, Species(x))))
    {}

    /// Construct an SpeciesListBase object with given another one.
    template<typename OtherData>
    SpeciesListBase(const SpeciesListBase<OtherData>& other)
     : m_species(other.data())
    {}

    /// Append a new species to the list of species.
    auto append(const Species& species)
    {
        m_species.push_back(species);
    }

    /// Return the internal collection of Species objects.
    auto data() const -> const Data&
    {
        return m_species;
    }

    /// Return the number of species in the collection.
    auto size() const -> Index
    {
        return m_species.size();
    }

    /// Return the Species object with given index.
    auto operator[](Index i) const -> const Species&
    {
        return m_species[i];
    }

    /// Return the index of the first species with given unique name or the number of species if not found.
    auto find(const String& name) const -> Index
    {
        return findWithName(name);
    }

    /// Return the index of the first species with given unique name or the number of species if not found.
    auto findWithName(const String& name) const -> Index
    {
        return indexfn(m_species, RKT_LAMBDA(s, s.name() == name));
    }

    /// Return the index of the first species with equivalent formula or the number of species if not found.
    auto findWithFormula(const ChemicalFormula& formula) const -> Index
    {
        return indexfn(m_species, RKT_LAMBDA(s, formula.equivalent(s.formula())));
    }

    /// Return the index of the first species with given substance name or the number of species if not found.
    auto findWithSubstance(const String& substance) const -> Index
    {
        return indexfn(m_species, RKT_LAMBDA(s, s.substance() == substance));
    }

    /// Return the index of the first species with given unique name or throw a runtime error if not found.
    auto index(const String& name) const -> Index
    {
        return indexWithName(name);
    }

    /// Return the index of the species with given unique name or throw a runtime error if not found.
    auto indexWithName(const String& name) const -> Index
    {
        const auto idx = findWithName(name);
        error(idx >= size(), "Could not find any Species object with name ", name, ".");
        return idx;
    }

    /// Return the index of the first species with equivalent formula or throw a runtime error if not found.
    auto indexWithFormula(const ChemicalFormula& formula) const -> Index
    {
        const auto idx = findWithFormula(formula);
        error(idx >= size(), "Could not find any Species object with formula ", formula.str(), ".");
        return idx;
    }

    /// Return the index of the first species with given substance name or throw a runtime error if not found.
    auto indexWithSubstance(const String& substance) const -> Index
    {
        const auto idx = findWithSubstance(substance);
        error(idx >= size(), "Could not find any Species object with substance ", substance, ".");
        return idx;
    }

    /// Return the species with given name.
    auto get(const String& name) const -> const Species&
    {
        return getWithName(name);
    }

    /// Return the species with given name.
    auto getWithName(const String& name) const -> const Species&
    {
        return m_species[indexWithName(name)];
    }

    /// Return the species with given formula.
    auto getWithFormula(const ChemicalFormula& formula) const -> const Species&
    {
        return m_species[indexWithFormula(formula)];
    }

    /// Return the species with given substance name.
    auto getWithSubstance(const String substance) const -> const Species&
    {
        return m_species[indexWithSubstance(substance)];
    }

    /// Return all species with given names.
    auto withNames(const StringList& names) const -> SpeciesList
    {
        return vectorize(names, RKT_LAMBDA(name, m_species[indexWithName(name)]));
    }

    /// Return all species with given substance formulas.
    auto withFormulas(const StringList& formulas) const -> SpeciesList
    {
        return vectorize(formulas, RKT_LAMBDA(formula, m_species[indexWithFormula(formula)]));
    }

    /// Return all species with given substance names.
    auto withSubstances(const StringList& substances) const -> SpeciesList
    {
        return vectorize(substances, RKT_LAMBDA(substance, m_species[indexWithSubstance(substance)]));
    }

    /// Return all species with given aggregate state.
    auto withAggregateState(AggregateState state) const -> SpeciesList
    {
        return filter(m_species, RKT_LAMBDA(s, s.aggregateState() == state));
    }

    /// Return all species with a given tag.
    auto withTag(String tag) const -> SpeciesList
    {
        if(tag.empty())
            return {};
        return filter(m_species, RKT_LAMBDA(s, contains(s.tags(), tag)));
    }

    /// Return all species without a given tag.
    auto withoutTag(String tag) const -> SpeciesList
    {
        if(tag.empty())
            return m_species;
        return filter(m_species, RKT_LAMBDA(s, !contains(s.tags(), tag)));
    }

    /// Return all species with given tags.
    auto withTags(const StringList& tags) const -> SpeciesList
    {
        if(tag.empty())
            return {};
        return filter(m_species, RKT_LAMBDA(s, contained(tags, s.tags())));
    }

    /// Return all species without given tags.
    auto withoutTags(const StringList& tags) const -> SpeciesList
    {
        if(tag.empty())
            return m_species;
        return filter(m_species, RKT_LAMBDA(s, !contained(tags, s.tags())));
    }

    /// Return all species with a certain elemental composition.
    /// This method filters the species composed of one or more given elements, as shown below.
    /// ~~~
    /// using namespace Reaktoro;
    /// SpeciesList specieslist("H2O H+ OH- H2 O2 Na+ Cl- NaCl CO2 HCO3- CO3-2 CH4");
    /// SpeciesList list1 = specieslist.withElements("H O");          // H2O H+ OH- H2 O2
    /// SpeciesList list2 = specieslist.withElements("H O C");        // H2O H+ OH- H2 O2 CO2 HCO3- CO3-2 CH4
    /// SpeciesList list4 = specieslist.withElements("H O Na Cl");    // H2O H+ OH- H2 O2 Na+ Cl- NaCl
    /// SpeciesList list5 = specieslist.withElements("H O Na Cl C");  // H2O H+ OH- H2 O2 Na+ Cl- NaCl CO2 HCO3- CO3-2 CH4
    /// ~~~
    /// @param symbols The element symbols of interest.
    /// @see SpeciesList::withElementsOf
    auto withElements(const StringList& symbols) const -> SpeciesList
    {
        return filter(m_species, RKT_LAMBDA(s, contained(s.elements().symbols(), symbols)));
    }

    /// Return all species with a certain elemental composition.
    /// This method extracts the elements from a given set of chemical formulas and
    /// then filters the species composed of one or more such elements, as shown below.
    /// ~~~
    /// using namespace Reaktoro;
    /// SpeciesList specieslist("H2O H+ OH- H2 O2 Na+ Cl- NaCl CO2 HCO3- CO3-2 CH4");
    /// SpeciesList list1 = specieslist.withElementsOf("H2O");         // H2O H+ OH- H2 O2
    /// SpeciesList list2 = specieslist.withElementsOf("H2O CO2");     // H2O H+ OH- H2 O2 CO2 HCO3- CO3-2 CH4
    /// SpeciesList list3 = specieslist.withElementsOf("H2O NaCl");    // H2O H+ OH- H2 O2 Na+ Cl- NaCl
    /// SpeciesList list4 = specieslist.withElementsOf("H2O Na+ Cl-"); // H2O H+ OH- H2 O2 Na+ Cl- NaCl
    /// ~~~
    /// @param formulas The formulas of the species from which elements are extracted.
    /// @see SpeciesList::withElements
    auto withElementsOf(const StringList& formulas) const -> SpeciesList
    {
        Strings symbols;
        for(ChemicalFormula formula : formulas)
            symbols = merge(symbols, formula.symbols());
        return withElements(symbols);
    }

    /// Convert this SpeciesListBase object into its Data.
    operator Data() { return m_species; }

    /// Convert this SpeciesListBase object into its Data.
    operator Data() const { return m_species; }

    // Ensure other SpeciesListBase types are friend among themselves.
    template<typename DataOther>
    friend class SpeciesListBase;

private:
    /// The species stored in the list.
    Data m_species;

public:
    /// Construct an SpeciesList object with given begin and end iterators.
    template<typename InputIterator>
    SpeciesListBase(InputIterator begin, InputIterator end) : m_species(begin, end) {}

    /// Return begin const iterator of this SpeciesList instance (for STL compatibility reasons).
    auto begin() const { return data().begin(); }

    /// Return begin iterator of this SpeciesList instance (for STL compatibility reasons).
    auto begin() { return data().begin(); }

    /// Return end const iterator of this SpeciesList instance (for STL compatibility reasons).
    auto end() const { return data().end(); }

    /// Return end iterator of this SpeciesList instance (for STL compatibility reasons).
    auto end() { return data().end(); }

    /// Append a new Species at the back of the container (for STL compatibility reasons).
    auto push_back(const Species& species) -> void { append(species); }

    /// Insert a container of Species objects into this SpeciesList instance (for STL compatibility reasons).
    template<typename Iterator, typename InputIterator>
    auto insert(Iterator pos, InputIterator begin, InputIterator end) -> void { m_species.insert(pos, begin, end); }

    /// The type of the value stored in a SpeciesList (for STL compatibility reasons).
    using value_type = Species;
};

} // namespace Reaktoro
