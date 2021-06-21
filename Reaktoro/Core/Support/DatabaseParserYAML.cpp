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

#include "DatabaseParserYAML.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/ParseUtils.hpp>
#include <Reaktoro/Common/YAML.hpp>
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Models/ReactionThermoModelYAML.hpp>
#include <Reaktoro/Models/StandardThermoModelYAML.hpp>
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
        for(auto child : doc["Elements"])
            addElement(child);

        for(auto child : doc["Species"])
            addSpecies(child);
    }

    /// Return the yaml node among element nodes with given unique @p symbol.
    auto getElementNode(const String& symbol) -> yaml
    {
        const auto node = doc["Elements"];
        if(!node) return {};
        for(auto i = 0; i < node.size(); ++i)
            if(node[i].at("Symbol").as<String>() == symbol)
                return node[i];
        return {};
    }

    /// Return the yaml node among species nodes with given unique @p name.
    auto getSpeciesNode(const String& name) -> yaml
    {
        const auto node = doc["Species"];
        if(!node.IsDefined()) return {};
        for(auto i = 0; i < node.size(); ++i)
            if(node[i].at("Name").as<String>() == name)
                return node[i];
        return {};
    }

    /// Add a new element with given symbol. Check if an element yaml node exists with such symbol first.
    auto addElement(const String& symbol) -> Element
    {
        const auto node = getElementNode(symbol);
        if(!node.IsNull())
            return addElement(node); // create an element using info in the database file.
        Element element(symbol); // create an element using info from default elements in Elements.
        element_list.append(element);
        return element;
    }

    /// Add a new element with given yaml @p node.
    auto addElement(const yaml& node) -> Element
    {
        assert(!node.IsNull());
        Element::Attribs attribs;
        node.copyRequiredChildValueTo("Symbol", attribs.symbol);
        node.copyRequiredChildValueTo("MolarMass", attribs.molar_mass);
        node.copyOptionalChildValueTo("Name", attribs.name, {});
        node.copyOptionalChildValueTo("Tags", attribs.tags, {});
        Element element(attribs);
        element_list.append(element);
        return element;
    }

    /// Add a new species with given unique @p name. An yaml node for this species must exist.
    auto addSpecies(const String& name) -> Species
    {
        const auto node = getSpeciesNode(name);
        errorif(node.IsNull(), "Cannot create a Species object with "
            "name `", name, "` without an yaml node with its data. "
            "Are you sure this name is correct and there is a species with this name in the database?");
        return addSpecies(node);
    }

    /// Add a new species with given yaml @p node.
    auto addSpecies(const yaml& node) -> Species
    {
        assert(node.IsDefined());
        errorif(!node["Name"], "Missing `Name` specification in:\n\n", node.repr());
        errorif(!node["Formula"], "Missing `Formula` specification in:\n\n", node.repr());
        errorif(!node["AggregateState"], "Missing `AggregateState` specification in:\n\n", node.repr());
        errorif(!node["Elements"], "Missing `Elements` specification in:\n\n", node.repr());
        errorif(!node["FormationReaction"] && !node["StandardThermoModel"], "Missing `FormationReaction` or `StandardThermoModel` specification in:\n\n", node.repr());
        const String name = node["Name"];
        const auto idx = species_list.find(name);
        if(idx < species_list.size())
            return species_list[idx]; // Do not add a species that has already been added! Return existing one.
        Species::Attribs attribs;
        node.copyRequiredChildValueTo("Name", attribs.name);
        node.copyRequiredChildValueTo("Formula", attribs.formula);
        node.copyOptionalChildValueTo("Substance", attribs.substance, {});
        node.copyOptionalChildValueTo("Charge", attribs.charge, 0.0);
        node.copyRequiredChildValueTo("AggregateState", attribs.aggregate_state);
        errorif(attribs.aggregate_state == AggregateState::Undefined,
            "Unsupported AggregateState value `", node["AggregateState"].as<String>(), "` in:\n\n", node.repr(), "\n\n"
            "The supported values are given below:\n\n", supportedAggregateStateValues());
        attribs.elements = createElementalComposition(node);
        attribs.formation_reaction = createFormationReaction(node);
        attribs.std_thermo_model = createStandardThermoModel(node);
        attribs.tags = createTags(node);
        Species species(attribs);
        species_list.append(species);
        return species;
    }

    /// Create the elemental composition of the species in the given yaml @p node.
    auto createElementalComposition(const yaml& node) -> ElementalComposition
    {
        auto child = node["Elements"];
        assert(!child.IsNull());
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
            .withReactionThermoModel(createReactionThermoModel(child))
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
    auto createReactionThermoModel(const yaml& node) -> ReactionThermoModel
    {
        auto child = node["ReactionThermoModel"];
        errorif(!child, "Missing `ReactionThermoModel` specification in:\n\n", node.repr());
        return ReactionThermoModelYAML(child);
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

auto DatabaseParserYAML::elements() const -> ElementListConstRef
{
    return pimpl->element_list;
}

auto DatabaseParserYAML::species() const -> SpeciesListConstRef
{
    return pimpl->species_list;
}

DatabaseParserYAML::operator Database() const
{
    return Database(elements(), species());
}

} // namespace Reaktoro
