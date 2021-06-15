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

#include "DatabaseDecoderYAML.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/ParseUtils.hpp>
#include <Reaktoro/Common/YAML.hpp>
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Models/StandardThermoModelYAML.hpp>
#include <Reaktoro/Serialization/Common.YAML.hpp>
#include <Reaktoro/Serialization/Core.YAML.hpp>
#include <Reaktoro/Serialization/Models.YAML.hpp>

namespace Reaktoro {

struct DatabaseDecoderYAML::Impl
{
    /// The Species objects in the database.
    SpeciesList species_list;

    /// The Element objects in the database.
    ElementList element_list;

    /// The document in yaml format with the database contents.
    yaml doc;

    /// Construct a default DatabaseDecoderYAML::Impl object.
    Impl()
    {}

    /// Construct a DatabaseDecoderYAML::Impl object with given YAML node.
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

    /// Add a new element with given symbol. Check if an element yaml node exists with such symbol first.
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
        node.copyRequiredChildValueTo("Symbol", attribs.symbol);
        node.copyRequiredChildValueTo("MolarMass", attribs.molar_mass);
        node.copyOptionalChildValueTo("Name", attribs.name);
        node.copyOptionalChildValueTo("Tags", attribs.tags);
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
        node.copyRequiredChildValueTo("Name", attribs.name);
        node.copyRequiredChildValueTo("Formula", attribs.formula);
        node.copyOptionalChildValueTo("Substance", attribs.substance);
        node.copyOptionalChildValueTo("Charge", attribs.charge);
        node.copyOptionalChildValueTo("AggregateState", attribs.aggregate_state);
        createElementalComposition(node, attribs.elements);
        createFormationReaction(node, attribs.formation_reaction);
        createStandardThermoModel(node, attribs.std_thermo_model);
        node.copyOptionalChildValueTo("Tags", attribs.tags);
        Species species(attribs);
        species_list.append(species);
        return species;
    }

    /// Create the elemental composition of the species in the given yaml @p node.
    auto createElementalComposition(const yaml& node, Optional<ElementalComposition>& elements) -> void
    {
        auto child = node["Elements"];
        if(!child) return; // there is no explicit info about the elements of the species (this is in the chemical formula)
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
        elements = ElementalComposition(pairs);
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

DatabaseDecoderYAML::DatabaseDecoderYAML()
: pimpl(new Impl())
{}

DatabaseDecoderYAML::DatabaseDecoderYAML(const DatabaseDecoderYAML& other)
: pimpl(new Impl(*other.pimpl))
{}

DatabaseDecoderYAML::DatabaseDecoderYAML(const yaml& doc)
: pimpl(new Impl(doc))
{}

DatabaseDecoderYAML::~DatabaseDecoderYAML()
{}

auto DatabaseDecoderYAML::operator=(DatabaseDecoderYAML other) -> DatabaseDecoderYAML&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

DatabaseDecoderYAML::operator Database() const
{
    return Database(pimpl->element_list, pimpl->species_list);
}

} // namespace Reaktoro
