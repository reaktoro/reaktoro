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

#include "DatabaseParser.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/YAML.hpp>
#include <Reaktoro/Models/StandardThermoModelYAML.hpp>
#include <Reaktoro/Serialization/Common.yaml.hpp>
#include <Reaktoro/Serialization/Core.yaml.hpp>
#include <Reaktoro/Serialization/Models.yaml.hpp>

namespace Reaktoro {

struct DatabaseParser::Impl
{
    /// The Species objects in the database.
    SpeciesList species_list;

    /// The Element objects in the database.
    ElementList element_list;

    /// The document in yaml format with the database contents.
    yaml doc;

    /// Construct a default DatabaseParser::Impl object.
    Impl()
    {}

    /// Construct a DatabaseParser::Impl object with given YAML node.
    Impl(const yaml& doc)
    {

    }

    /// Return the yaml node among element nodes with given unique @p symbol.
    auto getElementNode(const String& symbol) -> yaml
    {
        const auto node = doc["Elements"];
        if(!node.IsDefined()) return node;
        for(auto i = 0; i < node.size(); ++i)
            if(node[i].at("Symbol").as<std::string>() == symbol)
                return node[i];
        return {};
    }

    /// Return the yaml node among species nodes with given unique @p name.
    auto getSpeciesNode(const String& name) -> yaml
    {
        const auto node = doc["Species"];
        if(!node.IsDefined()) return node;
        for(auto i = 0; i < node.size(); ++i)
            if(node[i].at("Name").as<std::string>() == name)
                return node[i];
        return {};
    }

    /// Add a new element with given symbol. Check if an element yaml node exists if such symbol first.
    auto addElement(const String& symbol) -> Element
    {
        const auto node = getElementNode(symbol);
        if(node.IsDefined())
            return addElement(node); // create an element using info in the database file.
        Element element(symbol); // create an element using info from periodic table.
        element_list.append(element);
        return element;
    }

    /// Add a new element with given yaml @p node.
    auto addElement(const yaml& node) -> Element
    {
        assert(node.IsDefined());
        Element::Attribs attribs;
        node.required("Symbol", attribs.symbol);
        node.required("Name", attribs.name);
        node.required("AtomicNumber", attribs.atomic_number);
        node.required("AtomicWeight", attribs.atomic_weight);
        node.required("Electronegativity", attribs.electronegativity);
        node.required("Tags", attribs.tags);
        Element element(attribs);
        element_list.append(element);
        return element;
    }

    /// Add a new species with given unique @p name. An yaml node for this species must exist.
    auto addSpecies(const String& name) -> Species
    {
        const auto node = getSpeciesNode(name);
        errorif(!node.IsDefined(), "Cannot create a Species object with name `", name, "` without an yaml node with its data.");
        return addSpecies(node);
    }

    /// Add a new species with given yaml @p node.
    auto addSpecies(const yaml& node) -> Species
    {
        assert(node.IsDefined());
        Species::Attribs attribs;
        node.required("Name", attribs.name);
        node.required("Formula", attribs.formula);
        node.optional("Substance", attribs.substance);
        node.optional("Charge", attribs.charge);
        node.optional("AggregateState", attribs.aggregate_state);
        createElementalComposition(node, attribs.elements);
        createFormationReaction(node, attribs.formation_reaction);
        createStandardThermoModel(node, attribs.std_thermo_model);
        Species species(attribs);
        species_list.append(species);
        return species;
    }

    /// Create the elemental composition of the species in the given yaml @p node.
    auto createElementalComposition(const yaml& node, Optional<ElementalComposition>& elements) -> void
    {
        if(!node) return; // there is no explicit info about the elements of the species (this is in the chemical formula)
        Map<Element, double> elementMap;
        const auto pairs = node.as<std::map<String, double>>();
        for(const auto& [symbol, coeff] : pairs)
        {
            const auto idx = element_list.find(symbol);
            if(idx < element_list.size())
                elementMap.emplace(element_list[idx], coeff);
            else
            {
                const auto new_element = addElement(symbol);
                elementMap.emplace(new_element, coeff);
            }
        }
        elements = ElementalComposition(elementMap);
    }

    /// Create a standard thermodynamic model with given yaml @p node for a species.
    auto createStandardThermoModel(const yaml& node, Optional<StandardThermoModel>& model) -> void
    {
        auto child = node["StandardThermoModel"];
        if(!child) return;
        model = StandardThermoModelYAML(child);
    }

    /// Create a formation reaction with given yaml @p node for a species.
    auto createFormationReaction(const yaml& node, Optional<FormationReaction>& reaction) -> void
    {
        auto child = node["FormationReaction"];
        if(!child) return;
        reaction = FormationReaction()
            .withReactants(createReactants(child))
            .withReactionThermoModel(createReactionThermoModel(child))
            // .withProductStandardVolumeModel(createStandardVolumeModel(child))
            ;
    }

    /// Create the list of reactant species in given formation reaction yaml @p node.
    auto createReactants(const yaml& node) -> Pairs<Species, double>
    {
        Pairs<Species, double> reactants;

        auto pairs = node.as<std::map<String, double>>();
        for(const auto& [name, coeff] : pairs)
        {
            const auto idx = species_list.find(name);
            if(idx < species_list.size())
                reactants.emplace_back(species_list[idx], coeff);
            else
            {
                const auto new_species = addSpecies(name);
                reactants.emplace_back(new_species, coeff);
            }
        }
        return reactants;
    }

    /// Create a standard thermodynamic model with given yaml @p node for a species.
    auto createReactionThermoModel(const yaml& child) -> ReactionThermoModel
    {
        return {};
    }

    // auto createStandardVolumeModel(const yaml& child) -> StandardVolumeModel
    // {

    // }

};

DatabaseParser::DatabaseParser()
: pimpl(new Impl())
{}

DatabaseParser::DatabaseParser(const DatabaseParser& other)
: pimpl(new Impl(*other.pimpl))
{}

DatabaseParser::DatabaseParser(const yaml& doc)
: pimpl(new Impl(doc))
{}

DatabaseParser::~DatabaseParser()
{}

auto DatabaseParser::operator=(DatabaseParser other) -> DatabaseParser&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

} // namespace Reaktoro

