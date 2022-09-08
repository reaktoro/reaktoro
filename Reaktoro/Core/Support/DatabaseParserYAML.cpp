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

#include "DatabaseParserYAML.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/ParseUtils.hpp>
#include <Reaktoro/Common/YAML.hpp>
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelYAML.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelYAML.hpp>
#include <Reaktoro/Serialization/Common.YAML.hpp>
#include <Reaktoro/Serialization/Core.YAML.hpp>
#include <Reaktoro/Serialization/Models.YAML.hpp>

namespace Reaktoro {

struct DatabaseParserYAML::Impl
{
    /// The Species objects in the database.
    SpeciesList species_list;

    /// The Element objects in the database.
    ElementList element_list;

    /// The document in yaml format with the database contents.
    yaml doc;

    /// Construct a default DatabaseParserYAML::Impl object.
    Impl()
    {}

    /// Construct a DatabaseParserYAML::Impl object with given YAML node.
    Impl(const yaml& doc)
    : doc(doc)
    {
        errorif(!doc.IsMap(), "Could not understand your YAML database file with content:\n", doc.repr(), "\n",
            "Repeating the error message here in case the above printed content is too long.\n",
            "Could not understand your YAML database file with the above content.\n",
            "Are you forgetting to add the list of chemical species inside a Species YAML map?\n",
            "Please check other Reaktoro's YAML databases to identify what is not conforming.");

        for(auto const& child : doc["Elements"])
            addElement(child.first.as<String>(), child.second);

        for(auto const& child : doc["Species"])
            addSpecies(child.first.as<String>(), child.second);
    }

    /// Return the yaml node among element nodes with given unique @p symbol.
    auto getElementDetails(String const& symbol) -> yaml
    {
        const auto node = doc["Elements"];
        if(!node.IsDefined()) return {};
        auto const& child = node[symbol];
        return child ? yaml(child) : yaml();
    }

    /// Return the yaml node among species objects with given unique @p name.
    auto getSpeciesDetails(const String& name) -> yaml
    {
        const auto node = doc["Species"];
        if(!node.IsDefined()) return {};
        auto const& child = node[name];
        return child ? yaml(child) : yaml();
    }

    /// Add a new element with given symbol. Check if an element yaml node exists with such symbol first.
    auto addElement(String const& symbol) -> Element
    {
        const auto node = getElementDetails(symbol);
        if(!node.IsNull())
            return addElement(node); // create an element using info in the database file.
        Element element(symbol); // create an element using info from default elements in Elements.
        element_list.append(element);
        return element;
    }

    /// Add a new element with given `symbol` and `details`.
    auto addElement(String const& symbol, const yaml& details) -> Element
    {
        errorif(!details.IsMap(), "Expecting the details of an element as an object, but got instead:\n\n", details.repr());
        Element::Attribs attribs;
        attribs.symbol = symbol;
        details.copyRequiredChildValueTo("MolarMass", attribs.molar_mass);
        details.copyOptionalChildValueTo("Name", attribs.name, {});
        details.copyOptionalChildValueTo("Tags", attribs.tags, {});
        Element element(attribs);
        element_list.append(element);
        return element;
    }

    /// Add a new species with given unique @p name. An yaml node for this species must exist.
    auto addSpecies(String const& name) -> Species
    {
        const auto node = getSpeciesDetails(name);
        errorif(node.IsNull(), "Cannot create a Species object with "
            "name `", name, "` without an yaml node with its data. "
            "Are you sure this name is correct and there is a species with this name in the database?");
        return addSpecies(name, node);
    }

    /// Add a new species with given `name` and `details`.
    auto addSpecies(String const& name, yaml const& details) -> Species
    {
        errorif(!details.IsMap(), "Expecting the details of a species as an object, but got instead:\n\n", details.repr());
        errorif(!details["Formula"], "Missing `Formula` specification in:\n\n", details.repr());
        errorif(!details["AggregateState"], "Missing `AggregateState` specification in:\n\n", details.repr());
        errorif(!details["Elements"], "Missing `Elements` specification in:\n\n", details.repr(), "\n",
            "Please assign `Elements: null` if this species does not have chemical elements (e.g., e-, which may be represented with only `Charge: -1`).");
        errorif(!details["FormationReaction"] && !details["StandardThermoModel"], "Missing `FormationReaction` or `StandardThermoModel` specification in:\n\n", details.repr());
        const auto idx = species_list.find(name);
        if(idx < species_list.size())
            return species_list[idx]; // Do not add a species that has already been added! Return existing one.
        Species::Attribs attribs;
        attribs.name = name;
        details.copyRequiredChildValueTo("Formula", attribs.formula);
        details.copyOptionalChildValueTo("Substance", attribs.substance, {});
        details.copyOptionalChildValueTo("Charge", attribs.charge, 0.0);
        details.copyRequiredChildValueTo("AggregateState", attribs.aggregate_state);
        errorif(attribs.aggregate_state == AggregateState::Undefined,
            "Unsupported AggregateState value `", details["AggregateState"].as<String>(), "` in:\n\n", details.repr(), "\n\n"
            "The supported values are given below:\n\n", supportedAggregateStateValues());
        attribs.elements = createElementalComposition(details);
        attribs.formation_reaction = createFormationReaction(details);
        attribs.std_thermo_model = createStandardThermoModel(details);
        attribs.tags = createTags(details);
        Species species(attribs);
        species_list.append(species);
        return species;
    }

    /// Create the elemental composition of the species in the given yaml @p node.
    auto createElementalComposition(const yaml& node) -> ElementalComposition
    {
        auto child = node["Elements"];
        if(child.IsNull()) return {}; // for example, e- may be described with only Charge: -1 and Elements: null
        Pairs<Element, double> pairs;
        const auto symbols_and_coeffs = parseNumberStringPairs(child.as<String>());
        for(const auto& [symbol, coeff] : symbols_and_coeffs)
        {
            const auto idx = element_list.find(symbol);
            if(idx < element_list.size())
                pairs.emplace_back(element_list[idx], coeff);
            else
            {
                const auto new_element = addElement(symbol);
                pairs.emplace_back(new_element, coeff);
            }
        }
        return ElementalComposition(pairs);
    }

    /// Create the vector of tags with given yaml @p node for a species.
    auto createTags(const yaml& node) -> Strings
    {
        auto child = node["Tags"];
        if(!child.IsDefined()) return {};
        if(child.IsSequence()) return child.as<Strings>();
        return split(child.as<String>());
    }

    /// Create a standard thermodynamic model with given yaml @p node for a species.
    auto createStandardThermoModel(const yaml& node) -> StandardThermoModel
    {
        auto child = node["StandardThermoModel"];
        if(!child.IsDefined()) return {};
        return StandardThermoModelYAML(child);
    }

    /// Create a formation reaction with given yaml @p node for a species.
    auto createFormationReaction(const yaml& node) -> FormationReaction
    {
        auto child = node["FormationReaction"];
        if(!child.IsDefined()) return {};
        return FormationReaction()
            .withReactants(createReactants(child))
            .withReactionStandardThermoModel(createReactionStandardThermoModel(child))
            // .withProductStandardVolumeModel(createStandardVolumeModel(child))
            ;
    }

    /// Create the list of reactant species in given formation reaction yaml @p node.
    auto createReactants(const yaml& node) -> Pairs<Species, double>
    {
        auto child = node["Reactants"];
        const auto names_and_coeffs = parseNumberStringPairs(child.as<String>());
        Pairs<Species, double> reactants;
        for(const auto& [name, coeff] : names_and_coeffs)
        {
            const auto idx = species_list.find(name);
            if(idx < species_list.size()) // check if there is a Species object with current name
                reactants.emplace_back(species_list[idx], coeff); // if so, add it to the list of reactants
            else // otherwise, add a new Species object
            {
                const auto new_species = addSpecies(name); // Note potential recursivity: addSpecies may also need to call createReactants! This is needed in case reactions are defined recursively.
                reactants.emplace_back(new_species, coeff);
            }
        }
        return reactants;
    }

    /// Create a standard thermodynamic model with given yaml @p node for a species.
    auto createReactionStandardThermoModel(const yaml& node) -> ReactionStandardThermoModel
    {
        auto child = node["ReactionStandardThermoModel"];
        errorif(!child, "Missing `ReactionStandardThermoModel` specification in:\n\n", node.repr());
        return ReactionStandardThermoModelYAML(child);
    }

    // auto createStandardVolumeModel(const yaml& child) -> StandardVolumeModel
    // {

    // }
};

DatabaseParserYAML::DatabaseParserYAML()
: pimpl(new Impl())
{}

DatabaseParserYAML::DatabaseParserYAML(const DatabaseParserYAML& other)
: pimpl(new Impl(*other.pimpl))
{}

DatabaseParserYAML::DatabaseParserYAML(const yaml& doc)
: pimpl(new Impl(doc))
{}

DatabaseParserYAML::~DatabaseParserYAML()
{}

auto DatabaseParserYAML::operator=(DatabaseParserYAML other) -> DatabaseParserYAML&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto DatabaseParserYAML::elements() const -> const ElementList&
{
    return pimpl->element_list;
}

auto DatabaseParserYAML::species() const -> const SpeciesList&
{
    return pimpl->species_list;
}

DatabaseParserYAML::operator Database() const
{
    return Database(elements(), species());
}

} // namespace Reaktoro
