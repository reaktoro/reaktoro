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

#include "Core.yaml.hpp"

// Reaktoro includes
#include <Reaktoro/Core/ChemicalFormula.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/FormationReaction.hpp>
#include <Reaktoro/Core/Param.hpp>
#include <Reaktoro/Core/Params.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Serialization/Common.yaml.hpp>

namespace Reaktoro {

//=====================================================================================================================

REAKTORO_YAML_ENCODE_DEFINE(AggregateState)
{
    std::stringstream ss;
    ss << obj;
    node = ss.str();
}

REAKTORO_YAML_DECODE_DEFINE(AggregateState)
{
    obj = parseAggregateState(node.as<std::string>());
}

//=====================================================================================================================

REAKTORO_YAML_ENCODE_DEFINE(ChemicalFormula)
{
    node = obj.str();
}

REAKTORO_YAML_DECODE_DEFINE(ChemicalFormula)
{
    obj = ChemicalFormula(node.as<std::string>());
}

//=====================================================================================================================

REAKTORO_YAML_ENCODE_DEFINE(ChemicalSystem)
{
    errorif(true, "Implement REAKTORO_YAML_ENCODE_DEFINE(ChemicalSystem)");
}

REAKTORO_YAML_DECODE_DEFINE(ChemicalSystem)
{
    errorif(true, "Implement REAKTORO_YAML_DECODE_DEFINE(ChemicalSystem)");
}

//=====================================================================================================================

REAKTORO_YAML_ENCODE_DEFINE(Database)
{
    // ElementList elementlist;
    // auto child = node["Elements"];
    // if(child)
    // {
    //     for(auto i = 0; i < child.size(); i++)
    //         elementlist.append(child[i].as<Element>());
    // }
    // obj.
}

// REAKTORO_YAML_DECODE_DEFINE(Database)
// {
//     ElementList element_list;
//     SpeciesList species_list;

//     yaml node_species_list = node["Species"];
//     yaml node_element_list = node["Elements"];

//     auto appendElement = [&](const String& symbol)
//     {
//         Element new_element(symbol);
//         element_list.append(new_element);
//         return new_element;
//     };

//     createFormationReaction(child, attribs.formation_reaction);
//     createStandardThermoModel(child, attribs.std_thermo_model);

//     auto appendSpecies = [&](const String& name)
//     {
//         auto child = node_species_list[name];
//         errorif(!child, "Cannot construct database with missing data for species with name `", name, "`.");
//         Species::Attribs attribs;
//         child.optional("Name", attribs.name);
//         child.required("Formula", attribs.formula);
//         child.optional("Substance", attribs.substance);
//         child.optional("Elements", attribs.elements);
//         child.optional("Charge", attribs.charge);
//         child.optional("AggregateState", attribs.aggregate_state);
//         createFormationReaction(child, attribs.formation_reaction);
//         createStandardThermoModel(child, attribs.std_thermo_model);
//         // node.optional("StandardThermoModel", attribs.std_thermo_model);
//         // node.optional("StandardThermoModel", attribs.std_thermo_model);
//         node.optional("Tags", attribs.tags);

//         Element new_element(symbol);
//         element_list.append(new_element);
//         return new_element;
//     };

//     auto getElement = [&](const String& symbol)
//     {
//         const auto idx = element_list.find(symbol);
//         if(idx < element_list.size())
//             return element_list[idx];
//         return appendElement(symbol);
//     };

//     auto getSpecies = [&](const String& name)
//     {
//         const auto idx = species_list.find(name);
//         if(idx < species_list.size())
//             return species_list[idx];
//         return appendSpecies(name);
//     };

//     {
//         auto child = node["Elements"];
//         if(child)
//         {
//             for(auto i = 0; i < child.size(); i++)
//                 element_list.append(child[i].as<Element>());
//         }
//     }

//     {
//         auto child = node["Species"];
//         if(child)
//         {
//             for(auto i = 0; i < child.size(); i++)
//                 species_list.append(child[i].as<Species>());
//         }
//     }


// }

//=====================================================================================================================

REAKTORO_YAML_ENCODE_DEFINE(Element)
{
    node["Symbol"]            = obj.symbol();
    node["Name"]              = obj.name();
    node["AtomicNumber"]      = obj.atomicNumber();
    node["AtomicWeight"]      = obj.atomicWeight();
    node["Electronegativity"] = obj.electronegativity();
    node["Tags"]              = obj.tags();
}

REAKTORO_YAML_DECODE_DEFINE(Element)
{
    Element::Attribs attribs;
    node.required("Symbol", attribs.symbol);
    node.required("Name", attribs.name);
    node.required("AtomicNumber", attribs.atomic_number);
    node.required("AtomicWeight", attribs.atomic_weight);
    node.required("Electronegativity", attribs.electronegativity);
    node.required("Tags", attribs.tags);
    obj = Element(attribs);
}

//=====================================================================================================================

REAKTORO_YAML_ENCODE_DEFINE(ElementalComposition)
{
    errorif(true, "Implement REAKTORO_YAML_ENCODE_DEFINE(ElementalComposition)");
}

REAKTORO_YAML_DECODE_DEFINE(ElementalComposition)
{
    errorif(true, "Implement REAKTORO_YAML_DECODE_DEFINE(ElementalComposition)");
}

//=====================================================================================================================

REAKTORO_YAML_ENCODE_DEFINE(FormationReaction)
{
    errorif(true, "Implement REAKTORO_YAML_ENCODE_DEFINE(FormationReaction)");
}

REAKTORO_YAML_DECODE_DEFINE(FormationReaction)
{
    errorif(true, "Implement REAKTORO_YAML_DECODE_DEFINE(FormationReaction)");
}

//=====================================================================================================================

REAKTORO_YAML_ENCODE_DEFINE(Param)
{
    node = obj.value();
}

REAKTORO_YAML_DECODE_DEFINE(Param)
{
    obj = node.as<double>();
}

//=====================================================================================================================

REAKTORO_YAML_ENCODE_DEFINE(Params)
{
    node = obj.data();
}

REAKTORO_YAML_DECODE_DEFINE(Params)
{
    auto values = node.as<Vec<double>>();
    obj = Params(values.begin(), values.end());
}

//=====================================================================================================================

REAKTORO_YAML_ENCODE_DEFINE(Phase)
{
    errorif(true, "Implement REAKTORO_YAML_ENCODE_DEFINE(Phase)");
}

REAKTORO_YAML_DECODE_DEFINE(Phase)
{
    errorif(true, "Implement REAKTORO_YAML_DECODE_DEFINE(Phase)");
}

//=====================================================================================================================

REAKTORO_YAML_ENCODE_DEFINE(Species)
{
    errorif(true, "Implement REAKTORO_YAML_ENCODE_DEFINE(Species)");
}

REAKTORO_YAML_DECODE_DEFINE(Species)
{
    Species::Attribs attribs;
    node.required("Name", attribs.name);
    node.required("Formula", attribs.formula);
    node.optional("Substance", attribs.substance);
    node.optional("Elements", attribs.elements);
    node.optional("Charge", attribs.charge);
    node.optional("AggregateState", attribs.aggregate_state);
    node.optional("FormationReaction", attribs.formation_reaction);
    // node.optional("StandardThermoModel", attribs.std_thermo_model);
    node.optional("Tags", attribs.tags);
    obj = Species(attribs);
}

auto decodeElementalComposition(const yaml& node, const ElementList& elementlist) -> ElementalComposition
{
    const auto parsed = node.as<std::map<String, double>>();
    Map<Element, double> elements;
    for(const auto& [key, val] : parsed)
        elements.emplace(elementlist.getWithSymbol(key), val);

    // Map<Element, double> elements;
    // for(auto it = node.begin(); it != node.end(); ++it)
    // {
    //     const auto key = it->first.as<std::string>();
    //     const auto val = it->second.as<double>();
    //     elements.emplace(elementlist.getWithSymbol(key), val);
    // }
    return ElementalComposition(elements);
}

// auto decodeSpecies(const yaml& node, const ElementList& elementlist) -> Species
// {
//     Species::Attribs attribs;
//     node.optional("Name", attribs.name);
//     node.required("Formula", attribs.formula);
//     node.optional("Substance", attribs.substance);
//     // node.optional("Elements", attribs.elements);
//     // attribs.elements = decodeElementList(node.optional("Elements", attribs.elements);
//     node.optional("Charge", attribs.charge);
//     node.optional("AggregateState", attribs.aggregate_state);
//     node.optional("FormationReaction", attribs.formation_reaction);
//     // node.optional("StandardThermoModel", attribs.std_thermo_model);
//     node.optional("Tags", attribs.tags);
// }

} // namespace YAML
