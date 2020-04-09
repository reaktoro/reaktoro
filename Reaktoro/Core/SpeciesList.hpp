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
#include <optional>
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/StringList.hpp>
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
    explicit SpeciesList(std::vector<Species> species);

    /// Construct an SpeciesList object with given species formulas.
    explicit SpeciesList(StringList formulas);

    /// Append a new species to the list of species.
    auto append(Species species) -> void;

    /// Return the internal collection of Species objects.
    auto data() const -> const std::vector<Species>&;

    /// Return the number of species in the collection.
    auto size() const -> std::size_t;

    /// Return the Species object with given index.
    auto operator[](Index index) const -> const Species&;

    /// Return the index of the first species with given unique symbol.
    /// If there is no species with given unique symbol, return number of species.
    auto indexWithSymbol(std::string symbol) const -> Index;

    /// Return the index of the first species with given name.
    /// If there is no species with given name, return number of species.
    auto indexWithName(std::string name) const -> Index;

    /// Return the index of the first species with equivalent substance formula.
    /// If there is no species with given substance formula, return number of species.
    auto indexWithFormula(std::string formula) const -> Index;

    /// Return all species with given unique symbols.
    auto withSymbols(const StringList& symbols) const -> SpeciesList;

    /// Return all species with given names.
    auto withNames(const StringList& names) const -> SpeciesList;

    /// Return all species with given substance formulas.
    auto withFormulas(const StringList& formulas) const -> SpeciesList;

    /// Return all species with a given tag.
    auto withTag(std::string tag) const -> SpeciesList;

    /// Return all species without a given tag.
    auto withoutTag(std::string tag) const -> SpeciesList;

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
    /// @see Species::withElementsOf
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
    /// @see Species::withElements
    auto withElementsOf(const StringList& formulas) const -> SpeciesList;

    /// Return begin const iterator of this SpeciesList instance (for STL compatibility reasons).
    inline auto begin() const { return data().begin(); }

    /// Return begin iterator of this SpeciesList instance (for STL compatibility reasons).
    inline auto begin() { return data().begin(); }

    /// Return end const iterator of this SpeciesList instance (for STL compatibility reasons).
    inline auto end() const { return data().end(); }

    /// Return end iterator of this SpeciesList instance (for STL compatibility reasons).
    inline auto end() { return data().end(); }

    /// Append a new Species at the back of the container (for STL compatibility reasons).
    inline auto push_back(const Species& species) -> void { append(species); }

    /// The type of the value stored in a SpeciesList (for STL compatibility reasons).
    using value_type = Species;

private:
    /// The species stored in the list.
    std::vector<Species> m_species;
};

} // namespace Reaktoro
