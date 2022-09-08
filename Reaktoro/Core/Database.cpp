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

#include "Database.hpp"

// C++ includes
#include <fstream>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/ParseUtils.hpp>
#include <Reaktoro/Core/Support/DatabaseParserYAML.hpp>

namespace Reaktoro {

struct Database::Impl
{
    /// The Species objects in the database.
    SpeciesList species;

    /// The Element objects in the database.
    ElementList elements;

    /// The additional data in the database whose type is known at runtime only.
    Any attached_data;

    /// The symbols of all elements already in the database.
    Set<String> element_symbols;

    /// The names of all species already in the database.
    Set<String> species_names;

    /// The species in the database grouped in terms of their aggregate state
    Map<AggregateState, SpeciesList> species_with_aggregate_state;

    /// Add an element in the database.
    auto addElement(const Element& element) -> void
    {
        if(element_symbols.find(element.symbol()) == element_symbols.end())
        {
            elements.append(element);
            element_symbols.insert(element.symbol());
        }
    }

    /// Add a species in the database.
    auto addSpecies(Species newspecies) -> void
    {
        // Ensure a unique name is used when storing this species!
        auto name = newspecies.name();
        auto unique_name = name;
        while(species_names.find(unique_name) != species_names.end())
            unique_name = unique_name + "!"; // keep adding symbol ! to the name such as H2O, H2O!, H2O!!, H2O!!! if many H2O are given

        // Replace name if not unique.
        if(name != unique_name) {
            newspecies = newspecies.withName(unique_name);
            warning(true, "Species should have unique names in Database, but species ", name, " "
                "violates this rule. The unique name ", unique_name, " has been assigned instead.");
        }

        // Replace aggregate state to Aqueous if undefined. If this default
        // aggregate state option is not appropriate for a given scenario,
        // ensure then that the correct aggregate state of the species has
        // already been set in the Species object using method
        // Species::withAggregateState. This permits species names such as H2O,
        // CO2, H2, O2 as well as ions H+, OH-, HCO3- to be considered as
        // aqueous species, without requiring explicit aggregate state suffix
        // aq as in H2O(aq), CO2(aq), H2(aq). For all other species, ensure an
        // aggregate state suffix is provided, such as H2O(l), H2O(g),
        // CaCO3(s), Fe3O4(s).
        if(newspecies.aggregateState() == AggregateState::Undefined)
            newspecies = newspecies.withAggregateState(AggregateState::Aqueous);

        // Append the new Species object in the species container
        species.append(newspecies);

        // Update the list of unique species names.
        species_names.insert(newspecies.name());

        // Update the container of elements with the Element objects in this new species.
        for(auto&& [element, coeff] : newspecies.elements())
            addElement(element);

        // Add the new species in the group of species with same aggregate state
        species_with_aggregate_state[newspecies.aggregateState()].push_back(newspecies);
    }

    /// Construct a reaction with given equation.
    auto reaction(const String& equation) const -> Reaction
    {
        const auto pairs = parseReactionEquation(equation);
        Pairs<Species, double> reactants;
        for(const auto [name, coeff] : pairs)
            reactants.emplace_back(species.get(name), coeff);
        return Reaction().withEquation(reactants);
    }
};

Database::Database()
: pimpl(new Impl())
{}

Database::Database(const Database& other)
: pimpl(new Impl(*other.pimpl))
{}

Database::Database(const Vec<Element>& elements, const Vec<Species>& species)
: Database()
{
    for(const auto& x : elements)
        addElement(x);
    for(const auto& x : species)
        addSpecies(x);
}

Database::Database(const Vec<Species>& species)
: Database()
{
    for(const auto& x : species)
        addSpecies(x);
}

Database::~Database()
{}

auto Database::operator=(Database other) -> Database&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto Database::clear() -> void
{
    *pimpl = Database::Impl();
}

auto Database::addElement(const Element& element) -> void
{
    pimpl->addElement(element);
}

auto Database::addSpecies(const Species& species) -> void
{
    pimpl->addSpecies(species);
}

auto Database::addSpecies(const Vec<Species>& species) -> void
{
    for(const auto& x : species)
        addSpecies(x);
}

auto Database::attachData(const Any& data) -> void
{
    pimpl->attached_data = data;
}

auto Database::extend(const Database& other) -> void
{
    for(const auto& element : other.elements())
        addElement(element);

    for(const auto& species : other.species())
        addSpecies(species);

    // TODO: Replace Any by Map<String, Any> so that it becomes easier/more intuitive to unify different attached data to Database objects.
    // pimpl->attached_data = ???;
}

auto Database::elements() const -> const ElementList&
{
    return pimpl->elements;
}

auto Database::species() const -> const SpeciesList&
{
    return pimpl->species;
}

auto Database::speciesWithAggregateState(AggregateState option) const -> SpeciesList
{
    auto it = pimpl->species_with_aggregate_state.find(option);
    if(it == pimpl->species_with_aggregate_state.end())
        return {};
    return it->second;
}

auto Database::element(const String& symbol) const -> const Element&
{
    return elements().getWithSymbol(symbol);
}

auto Database::species(const String& name) const -> const Species&
{
    return species().getWithName(name);
}

auto Database::reaction(const String& equation) const -> Reaction
{
    return pimpl->reaction(equation);
}

auto Database::attachedData() const -> const Any&
{
    return pimpl->attached_data;
}

auto Database::fromFile(const String& path) -> Database
{
    std::ifstream file(path);
    errorif(!file.is_open(),
        "Could not open file `", path, "`. Ensure the given file path "
        "is relative to the directory where your application is RUNNING "
        "(not necessarily where the executable is located!). Alternatively, "
        "try a full path to the file (e.g., "
        "in Windows, `C:\\User\\username\\mydata\\mydatabase.yaml`, "
        "in Linux, `/home/username/mydata/mydatabase.yaml`).");
    return fromStream(file);
}

auto Database::fromContents(const String& contents) -> Database
{
    auto doc = yaml::parse(contents);
    DatabaseParserYAML dbparser(doc);
    return dbparser;
}

auto Database::fromStream(std::istream& stream) -> Database
{
    auto doc = yaml::parse(stream);
    DatabaseParserYAML dbparser(doc);
    return dbparser;
}

} // namespace Reaktoro

