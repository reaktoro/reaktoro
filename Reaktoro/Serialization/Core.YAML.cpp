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

#include "Core.YAML.hpp"

// Reaktoro includes
#include <Reaktoro/Core/ChemicalFormula.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/ElementList.hpp>
#include <Reaktoro/Core/FormationReaction.hpp>
#include <Reaktoro/Core/Param.hpp>
#include <Reaktoro/Core/Params.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>
#include <Reaktoro/Core/Support/DatabaseParserYAML.hpp>
#include <Reaktoro/Models/ReactionThermoModelYAML.hpp>
#include <Reaktoro/Models/StandardThermoModelYAML.hpp>
#include <Reaktoro/Serialization/Common.YAML.hpp>

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
    auto elements_node = node["Elements"];
    auto species_node = node["Species"];
    elements_node = obj.elements();
    species_node = obj.species();
}

REAKTORO_YAML_DECODE_DEFINE(Database)
{
    obj = DatabaseParserYAML(node);
    obj = DatabaseParserYAML(node);
}

//=====================================================================================================================

REAKTORO_YAML_ENCODE_DEFINE(Element)
{
    node["Symbol"]    = obj.symbol();
    node["MolarMass"] = obj.molarMass();
    node.appendIfNotDefault("Name", obj.name(), "");
    node.appendIfNotDefault("Tags", obj.tags(), {});
}

REAKTORO_YAML_DECODE_DEFINE(Element)
{
    Element::Attribs attribs;
    node.copyRequiredChildValueTo("Symbol", attribs.symbol);
    node.copyRequiredChildValueTo("MolarMass", attribs.molar_mass);
    node.copyOptionalChildValueTo("Name", attribs.name, {});
    node.copyOptionalChildValueTo("Tags", attribs.tags, {});
    obj = Element(attribs);
}

//=====================================================================================================================

REAKTORO_YAML_ENCODE_DEFINE(ElementList)
{
    for(const auto& element : obj)
        node.push_back(element);
}

REAKTORO_YAML_DECODE_DEFINE(ElementList)
{
    obj = node.as<Vec<Element>>();
}

//=====================================================================================================================

REAKTORO_YAML_ENCODE_DEFINE(ElementalComposition)
{
    node = obj.repr();
}

REAKTORO_YAML_DECODE_DEFINE(ElementalComposition)
{
    errorif(true, "Converting YAML to ElementalComposition is not supported directly."); // because only element symbols are present in YAML representation of ElementalComposition
}

//=====================================================================================================================

REAKTORO_YAML_ENCODE_DEFINE(FormationReaction)
{
    std::stringstream repr;
    auto i = 0;
    for(const auto& [species, coeff] : obj.reactants())
        repr << (i++ == 0 ? "" : " ") << coeff << ":" << species.name();

    node["Reactants"] = repr.str();
    if(obj.reactionThermoModel().initialized())
        node["ReactionThermoModel"] = obj.reactionThermoModel().serialize();
    if(obj.productStandardVolumeModel().initialized())
        node["StandardVolumeModel"] = obj.productStandardVolumeModel().serialize();
}

REAKTORO_YAML_DECODE_DEFINE(FormationReaction)
{
    errorif(true, "Converting YAML to FormationReaction is not supported directly."); // because only species names are present in YAML representation of Reactants
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

// REAKTORO_YAML_ENCODE_DEFINE(Params)
// {

// }

// REAKTORO_YAML_DECODE_DEFINE(Params)
// {

// }

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

REAKTORO_YAML_ENCODE_DEFINE(ReactionThermoModel)
{
    node = obj.serialize();
}

REAKTORO_YAML_DECODE_DEFINE(ReactionThermoModel)
{
    obj = ReactionThermoModelYAML(node);
}

//=====================================================================================================================

REAKTORO_YAML_ENCODE_DEFINE(Species)
{
    node["Name"] = obj.name();
    node["Formula"] = obj.formula();
    node["Substance"] = obj.substance();
    node["Elements"] = obj.elements();
    node.appendIfNotDefault("Charge", obj.charge(), 0.0);
    node["AggregateState"] = obj.aggregateState();
    if(obj.reaction().reactants().size())
        node["FormationReaction"] = obj.reaction();
    else node["StandardThermoModel"] = obj.standardThermoModel();
    node.appendIfNotDefault("Tags", obj.tags(), Strings{});
}

REAKTORO_YAML_DECODE_DEFINE(Species)
{
    errorif(true, "Converting YAML to Species is not supported directly."); // because only element symbols are present in YAML representation of Species
}

//=====================================================================================================================

REAKTORO_YAML_ENCODE_DEFINE(SpeciesList)
{
    for(const auto& species : obj)
        node.push_back(species);
}

REAKTORO_YAML_DECODE_DEFINE(SpeciesList)
{
    errorif(true, "Converting YAML to SpeciesList is not supported directly."); // because only element symbols are present in YAML representation of Species
}

//=====================================================================================================================

REAKTORO_YAML_ENCODE_DEFINE(StandardThermoModel)
{
    node = obj.serialize();
}

REAKTORO_YAML_DECODE_DEFINE(StandardThermoModel)
{
    obj = StandardThermoModelYAML(node);
}

//=====================================================================================================================

} // namespace YAML
