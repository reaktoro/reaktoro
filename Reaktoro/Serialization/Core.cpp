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

#include "Core.hpp"

// Reaktoro includes
#include <Reaktoro/Core/ChemicalFormula.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/ElementList.hpp>
#include <Reaktoro/Core/FormationReaction.hpp>
#include <Reaktoro/Core/Param.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>
#include <Reaktoro/Core/Support/DatabaseParser.hpp>
#include <Reaktoro/Models/StandardThermoModels.hpp>
#include <Reaktoro/Serialization/Common.hpp>
#include <Reaktoro/Serialization/Models/ReactionRateModels.hpp>
#include <Reaktoro/Serialization/Models/StandardThermoModels.hpp>

namespace Reaktoro {

//=====================================================================================================================

REAKTORO_DATA_ENCODE_DEFINE(ChemicalFormula)
{
    data = obj.str();
}

REAKTORO_DATA_DECODE_DEFINE(ChemicalFormula)
{
    obj = ChemicalFormula(data.asString());
}

//=====================================================================================================================

REAKTORO_DATA_ENCODE_DEFINE(ChemicalSystem)
{
    errorif(true, "Implement REAKTORO_DATA_ENCODE_DEFINE(ChemicalSystem)");
}

REAKTORO_DATA_DECODE_DEFINE(ChemicalSystem)
{
    errorif(true, "Implement REAKTORO_DATA_DECODE_DEFINE(ChemicalSystem)");
}

//=====================================================================================================================

REAKTORO_DATA_ENCODE_DEFINE(Database)
{
    auto elements_data = data.at("Elements");
    auto species_data = data.at("Species");

    for(auto const& element : obj.elements())
        elements_data[element.symbol()] = element;

    for(auto const& species : obj.species())
        species_data[species.name()] = species;
}

REAKTORO_DATA_DECODE_DEFINE(Database)
{
    obj = DatabaseParser(data);
}

//=====================================================================================================================

REAKTORO_DATA_ENCODE_DEFINE(Element)
{
    data["Symbol"]    = obj.symbol();
    data["MolarMass"] = obj.molarMass();
    if(obj.name().size()) data["Name"] = obj.name();
    if(obj.tags().size()) data["Tags"] = obj.tags();
}

REAKTORO_DATA_DECODE_DEFINE(Element)
{
    Element::Attribs attribs;
    data.at("Symbol").to(attribs.symbol);
    data.at("MolarMass").to(attribs.molar_mass);
    if(data.exists("Name")) data.at("Name").to(attribs.name);
    if(data.exists("Tags")) data.at("Tags").to(attribs.tags);
    obj = Element(attribs);
}

//=====================================================================================================================

REAKTORO_DATA_ENCODE_DEFINE(ElementList)
{
    for(auto const& element : obj)
        data.add(element);
}

REAKTORO_DATA_DECODE_DEFINE(ElementList)
{
    obj = data.as<Vec<Element>>();
}

//=====================================================================================================================

REAKTORO_DATA_ENCODE_DEFINE(ElementalComposition)
{
    data = obj.repr();
}

REAKTORO_DATA_DECODE_DEFINE(ElementalComposition)
{
    errorif(true, "Converting YAML to ElementalComposition is not supported directly."); // because only element symbols are present in YAML representation of ElementalComposition
}

//=====================================================================================================================

REAKTORO_DATA_ENCODE_DEFINE(FormationReaction)
{
    std::stringstream repr;
    auto i = 0;
    for(const auto& [species, coeff] : obj.reactants())
        repr << (i++ == 0 ? "" : " ") << coeff << ":" << species.name();

    data["Reactants"] = repr.str();
    if(obj.reactionThermoModel().initialized())
        data["ReactionStandardThermoModel"] = obj.reactionThermoModel().serialize();
    if(obj.productStandardVolumeModel().initialized())
        data["StandardVolumeModel"] = obj.productStandardVolumeModel().serialize();
}

REAKTORO_DATA_DECODE_DEFINE(FormationReaction)
{
    errorif(true, "Converting YAML to FormationReaction is not supported directly."); // because only species names are present in YAML representation of Reactants
}

//=====================================================================================================================

REAKTORO_DATA_ENCODE_DEFINE(Phase)
{
    errorif(true, "Implement REAKTORO_DATA_ENCODE_DEFINE(Phase)");
}

REAKTORO_DATA_DECODE_DEFINE(Phase)
{
    errorif(true, "Implement REAKTORO_DATA_DECODE_DEFINE(Phase)");
}

//=====================================================================================================================

REAKTORO_DATA_ENCODE_DEFINE(ReactionStandardThermoModel)
{
    data = obj.serialize();
}

REAKTORO_DATA_DECODE_DEFINE(ReactionStandardThermoModel)
{
    auto createModel = [](Data const& data) -> ReactionStandardThermoModel
    {
        errorif(!data.isDict(),
            "Expecting a dictionary containing a single key-value pair in the Data object:\n", data.repr());

        const auto model = data.asDict();

        errorif(model.size() != 1,
            "Expecting only one key-value pair in the Data object:\n", data.repr());

        auto const& name = model.front().first;
        auto const& params = model.front().second;

        if(name == "ConstLgK")
            return ReactionStandardThermoModelConstLgK(params.as<ReactionStandardThermoModelParamsConstLgK>());
        if(name == "GemsLgK")
            return ReactionStandardThermoModelGemsLgK(params.as<ReactionStandardThermoModelParamsGemsLgK>());
        if(name == "PhreeqcLgK")
            return ReactionStandardThermoModelPhreeqcLgK(params.as<ReactionStandardThermoModelParamsPhreeqcLgK>());
        if(name == "VantHoff")
            return ReactionStandardThermoModelVantHoff(params.as<ReactionStandardThermoModelParamsVantHoff>());

        errorif(true, "Cannot create a ReactionStandardThermoModel object with "
            "unsupported model name `", name, "` in Data object:\n", data.repr());

        return {};
    };

    obj = createModel(data);
}

//=====================================================================================================================

REAKTORO_DATA_ENCODE_DEFINE(Species)
{
    data["Name"] = obj.name();
    data["Formula"] = obj.formula();
    data["Substance"] = obj.substance();
    data["Elements"] = obj.elements();
    if(obj.charge() != 0.0) data["Charge"] = obj.charge();
    data["AggregateState"] = obj.aggregateState();
    if(obj.reaction().reactants().size())
        data["FormationReaction"] = obj.reaction();
    else data["StandardThermoModel"] = obj.standardThermoModel();
    if(obj.tags().size()) data["Tags"] = obj.tags();
}

REAKTORO_DATA_DECODE_DEFINE(Species)
{
    errorif(true, "Converting YAML to Species is not supported directly."); // because only element symbols are present in YAML representation of Species
}

//=====================================================================================================================

REAKTORO_DATA_ENCODE_DEFINE(SpeciesList)
{
    for(auto const& species : obj)
        data.add(species);
}

REAKTORO_DATA_DECODE_DEFINE(SpeciesList)
{
    errorif(true, "Converting YAML to SpeciesList is not supported directly."); // because only element symbols are present in YAML representation of Species
}

//=====================================================================================================================

REAKTORO_DATA_ENCODE_DEFINE(StandardThermoModel)
{
    data = obj.serialize();
}

REAKTORO_DATA_DECODE_DEFINE(StandardThermoModel)
{
    auto createModel = [](Data const& data) -> StandardThermoModel
    {
        errorif(!data.isDict(),
            "Expecting a dictionary containing a single key-value pair in the Data object:\n", data.repr());

        const auto model = data.asDict();

        errorif(model.size() != 1,
            "Expecting only one key-value pair in the Data object:\n", data.repr());

        auto const& name = model.front().first;
        auto const& params = model.front().second;

        if(name == "Constant")
            return StandardThermoModelConstant(params.as<StandardThermoModelParamsConstant>());
        if(name == "HKF")
            return StandardThermoModelHKF(params.as<StandardThermoModelParamsHKF>());
        if(name == "HollandPowell")
            return StandardThermoModelHollandPowell(params.as<StandardThermoModelParamsHollandPowell>());
        if(name == "Interpolation")
            return StandardThermoModelInterpolation(params.as<StandardThermoModelParamsInterpolation>());
        if(name == "MaierKelley")
            return StandardThermoModelMaierKelley(params.as<StandardThermoModelParamsMaierKelley>());
        if(name == "MineralHKF")
            return StandardThermoModelMineralHKF(params.as<StandardThermoModelParamsMineralHKF>());
        if(name == "WaterHKF")
            return StandardThermoModelWaterHKF(params.as<StandardThermoModelParamsWaterHKF>());
        if(name == "Nasa")
            return StandardThermoModelNasa(params.as<StandardThermoModelParamsNasa>());

        errorif(true, "Cannot create a StandardThermoModel object with "
            "unsupported model name `", name, "` in Data object:\n", data.repr());

        return {};
    };

    obj = createModel(data);
}

//=====================================================================================================================

} // namespace Reaktoro
