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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/ElementList.hpp>
#include <Reaktoro/Core/Species.hpp>

namespace Reaktoro {

/// A type used as a collection of species.
class SpeciesList
{
public:
    /// Construct a default SpeciesList object.
    SpeciesList();

    /// Construct an SpeciesList object with given species.
    SpeciesList(std::initializer_list<Species> species);

    /// Construct an SpeciesList object with given species.
    SpeciesList(const Vec<Species>& species);

    /// Construct an SpeciesList object with given species formulas.
    SpeciesList(const StringList& formulas);

    /// Append a new species to the list of species.
    auto append(const Species& species) -> void;

    /// Return the internal collection of Species objects.
    auto data() const -> const Vec<Species>&;

    /// Return true if there are no species in the collection.
    auto empty() const -> bool;

    /// Return the number of species in the collection.
    auto size() const -> Index;

    /// Return the elements that compose the species in the collection sorted in ascending order of molar mass.
    auto elements() const -> ElementList;

    /// Return the Species object with given index.
    auto operator[](Index i) const -> const Species&;

    /// Return the Species object with given index.
    auto operator[](Index i) -> Species&;

    /// Return the index of the first species with given unique name or the number of species if not found.
    auto find(const String& name) const -> Index;

    /// Return the index of the first species with given unique name or the number of species if not found.
    auto findWithName(const String& name) const -> Index;

    /// Return the index of the first species with equivalent formula or the number of species if not found.
    auto findWithFormula(const ChemicalFormula& formula) const -> Index;

    /// Return the index of the first species with given substance name or the number of species if not found.
    auto findWithSubstance(const String& substance) const -> Index;

    /// Return the index of the first species with given unique name or throw a runtime error if not found.
    auto index(const String& name) const -> Index;

    /// Return the index of the species with given unique name or throw a runtime error if not found.
    auto indexWithName(const String& name) const -> Index;

    /// Return the index of the first species with equivalent formula or throw a runtime error if not found.
    auto indexWithFormula(const ChemicalFormula& formula) const -> Index;

    /// Return the index of the first species with given substance name or throw a runtime error if not found.
    auto indexWithSubstance(const String& substance) const -> Index;

    /// Return the species with a given name.
    auto get(const String& name) const -> const Species&;

    /// Return the species with a given name.
    auto getWithName(const String& name) const -> const Species&;

    /// Return the species with a given formula.
    auto getWithFormula(const ChemicalFormula& formula) const -> const Species&;

    /// Return the species with a given substance name.
    auto getWithSubstance(const String substance) const -> const Species&;

    /// Return all species with given names.
    auto withNames(const StringList& names) const -> SpeciesList;

    /// Return all species with given substance formulas.
    auto withFormulas(const StringList& formulas) const -> SpeciesList;

    /// Return all species with given substance names.
    auto withSubstances(const StringList& substances) const -> SpeciesList;

    /// Return all species with a given aggregate state.
    auto withAggregateState(AggregateState state) const -> SpeciesList;

    /// Return all species with a given charge.
    auto withCharge(real value) const -> SpeciesList;

    /// Return all species with a given tag.
    auto withTag(String tag) const -> SpeciesList;

    /// Return all species without a given tag.
    auto withoutTag(String tag) const -> SpeciesList;

    /// Return all species with given tags.
    auto withTags(const StringList& tags) const -> SpeciesList;

    /// Return all species without given tags.
    auto withoutTags(const StringList& tags) const -> SpeciesList;

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
    auto withElements(const StringList& symbols) const -> SpeciesList;

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
    auto withElementsOf(const StringList& formulas) const -> SpeciesList;

    /// Convert this SpeciesList object into its Data.
    operator Vec<Species>&();

    /// Convert this SpeciesList object into its Data.
    operator Vec<Species>const&() const;

private:
    /// The species stored in the list.
    Vec<Species> m_species;

public:
    /// Construct an SpeciesList object with given begin and end iterators.
    template<typename InputIterator>
    SpeciesList(InputIterator begin, InputIterator end) : m_species(begin, end) {}

    /// Return begin const iterator of this SpeciesList instance (for STL compatibility reasons).
    auto begin() const { return m_species.begin(); }

    /// Return begin iterator of this SpeciesList instance (for STL compatibility reasons).
    auto begin() { return m_species.begin(); }

    /// Return end const iterator of this SpeciesList instance (for STL compatibility reasons).
    auto end() const { return m_species.end(); }

    /// Return end iterator of this SpeciesList instance (for STL compatibility reasons).
    auto end() { return m_species.end(); }

    /// Append a new Species at the back of the container (for STL compatibility reasons).
    auto push_back(const Species& species) -> void { append(species); }

    /// Insert a container of Species objects into this SpeciesList instance (for STL compatibility reasons).
    template<typename Iterator, typename InputIterator>
    auto insert(Iterator pos, InputIterator begin, InputIterator end) -> void { m_species.insert(pos, begin, end); }

    /// The type of the value stored in a SpeciesList (for STL compatibility reasons).
    using value_type = Species;
};

/// Return the concatenation of two SpeciesList objects.
auto operator+(const SpeciesList& a, const SpeciesList& b) -> SpeciesList;

} // namespace Reaktoro
