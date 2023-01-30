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

#include "ActivityModels.hpp"

// Reaktoro includes
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Models/ActivityModels.hpp>

namespace Reaktoro {

//======================================================================
// ActivityModelParams Types
//======================================================================

using PitzerCorrectionModel         = ActivityModelParamsPitzer::CorrectionModel;
using PitzerInteractionParamAttribs = ActivityModelParamsPitzer::InteractionParamAttribs;
using PitzerAlphaParamAttribs       = ActivityModelParamsPitzer::AlphaParamAttribs;

auto pitzerCorrectionModelToString(PitzerCorrectionModel const& modeltype) -> String
{
    switch(modeltype)
    {
        case PitzerCorrectionModel::Phreeqc:            return "Phreeqc";
        case PitzerCorrectionModel::HeMorse1993:        return "HeMorse1993";
        case PitzerCorrectionModel::Dai2013:            return "Dai2013";
        case PitzerCorrectionModel::Dai2014:            return "Dai2014";
        case PitzerCorrectionModel::ChristovMoller2004: return "ChristovMoller2004";
        case PitzerCorrectionModel::Holmes1987:         return "Holmes1987";
        case PitzerCorrectionModel::Pitzer1984:         return "Pitzer1984";
        case PitzerCorrectionModel::PalabanPitzer1987:  return "PalabanPitzer1987";
        case PitzerCorrectionModel::Polya2001:          return "Polya2001";
        case PitzerCorrectionModel::LiDuan2007:         return "LiDuan2007";
    }
    errorif(true, "Unable to encode Pitzer parameter (T, P) correction model `", modeltype, "` to a string.");
    return "";
};

auto pitzerCorrectionModelFromString(String const& modeltype) -> PitzerCorrectionModel
{
    if(modeltype == "Phreeqc")             return PitzerCorrectionModel::Phreeqc;
    if(modeltype == "HeMorse1993")         return PitzerCorrectionModel::HeMorse1993;
    if(modeltype == "Dai2013")             return PitzerCorrectionModel::Dai2013;
    if(modeltype == "Dai2014")             return PitzerCorrectionModel::Dai2014;
    if(modeltype == "ChristovMoller2004")  return PitzerCorrectionModel::ChristovMoller2004;
    if(modeltype == "Holmes1987")          return PitzerCorrectionModel::Holmes1987;
    if(modeltype == "Pitzer1984")          return PitzerCorrectionModel::Pitzer1984;
    if(modeltype == "PalabanPitzer1987")   return PitzerCorrectionModel::PalabanPitzer1987;
    if(modeltype == "Polya2001")           return PitzerCorrectionModel::Polya2001;
    if(modeltype == "LiDuan2007")          return PitzerCorrectionModel::LiDuan2007;

    errorif(true, "Unable to decode Pitzer parameter (T, P) correction model from string `", modeltype, "`.");
    return {};
};

REAKTORO_DATA_ENCODE_DEFINE(ActivityModelParamsPitzer)
{
    auto encodeInteractionParamAttribs = [](PitzerInteractionParamAttribs const& attribs) -> Data
    {
        Data data;
        data.add(pitzerCorrectionModelToString(attribs.model));
        for(auto const& formula : attribs.formulas)
            data.add(formula);
        for(auto const& coeff : attribs.coefficients)
            data.add(coeff);
        return data;
    };

    auto encodeAlphaParamAttribs = [](PitzerAlphaParamAttribs const& attribs) -> Data
    {
        Data data;
        data.add(attribs.formula1);
        data.add(attribs.formula2);
        data.add(attribs.value);
        return data;
    };

    for(auto const& param : obj.beta0)
        data["beta0"].add(encodeInteractionParamAttribs(param));

    for(auto const& param : obj.beta1)
        data["beta1"].add(encodeInteractionParamAttribs(param));

    for(auto const& param : obj.beta2)
        data["beta2"].add(encodeInteractionParamAttribs(param));

    for(auto const& param : obj.Cphi)
        data["Cphi"].add(encodeInteractionParamAttribs(param));

    for(auto const& param : obj.theta)
        data["theta"].add(encodeInteractionParamAttribs(param));

    for(auto const& param : obj.psi)
        data["psi"].add(encodeInteractionParamAttribs(param));

    for(auto const& param : obj.lambda)
        data["lambda"].add(encodeInteractionParamAttribs(param));

    for(auto const& param : obj.zeta)
        data["zeta"].add(encodeInteractionParamAttribs(param));

    for(auto const& param : obj.mu)
        data["mu"].add(encodeInteractionParamAttribs(param));

    for(auto const& param : obj.eta)
        data["eta"].add(encodeInteractionParamAttribs(param));

    for(auto const& param : obj.alpha1)
        data["alpha1"].add(encodeAlphaParamAttribs(param));

    for(auto const& param : obj.alpha2)
        data["alpha2"].add(encodeAlphaParamAttribs(param));
}

REAKTORO_DATA_DECODE_DEFINE(ActivityModelParamsPitzer)
{
    errorifnot(data.isDict(), "While decoding the data shown next, I was expecting a YAML map or a JSON object but got instead:\n", data.dumpYaml());

    auto extractCorrectionModel = [](Data const& row, Vec<Data> const& values) -> PitzerCorrectionModel
    {
        errorifnot(values.size() > 0 && values[0].isString(), "Expecting the 1st entry in the following list to be a supported representative string for a ActivityModelParamsPitzer::CorrectionModel value:\n", row.dumpYaml());
        return pitzerCorrectionModelFromString(values[0].asString());
    };

    auto extractChemicalFormulasBinary = [](Data const& row, Vec<Data> const& values) -> Vec<ChemicalFormula>
    {
        errorifnot(values.size() > 1 && values[1].isString(), "Expecting the 2nd entry in the following list to be a valid chemical formula:\n", row.dumpYaml());
        errorifnot(values.size() > 2 && values[2].isString(), "Expecting the 3rd entry in the following list to be a valid chemical formula:\n", row.dumpYaml());
        return {values[1].asString(), values[2].asString()};
    };

    auto extractChemicalFormulasTernary = [](Data const& row, Vec<Data> const& values) -> Vec<ChemicalFormula>
    {
        errorifnot(values.size() > 1 && values[1].isString(), "Expecting the 2nd entry in the following list to be a valid chemical formula:\n", row.dumpYaml());
        errorifnot(values.size() > 2 && values[2].isString(), "Expecting the 3rd entry in the following list to be a valid chemical formula:\n", row.dumpYaml());
        errorifnot(values.size() > 3 && values[3].isString(), "Expecting the 4th entry in the following list to be a valid chemical formula:\n", row.dumpYaml());
        return {values[1].asString(), values[2].asString(), values[3].asString()};
    };

    auto extractCoefficientsBinary = [](Data const& row, Vec<Data> const& values) -> Vec<Param>
    {
        errorifnot(values.size() > 3, "Expecting at least one coefficient value in the following list:\n", row.dumpYaml());
        Vec<Param> coefficients;
        for(auto i = 3; i <= values.size(); ++i)
            coefficients.push_back(values[i].asParam());
        return coefficients;
    };

    auto extractCoefficientsTernary = [](Data const& row, Vec<Data> const& values) -> Vec<Param>
    {
        errorifnot(values.size() > 4, "Expecting at least one coefficient value in the following list:\n", row.dumpYaml());
        Vec<Param> coefficients;
        for(auto i = 4; i <= values.size(); ++i)
            coefficients.push_back(values[i].asParam());
        return coefficients;
    };

    auto decodeInteractionParamAttribsBinary = [&](Data const& row) -> PitzerInteractionParamAttribs
    {
        errorifnot(row.isList(), "While decoding the data shown next, I was expecting a YAML/JSON list but got instead:\n", row.dumpYaml());
        const auto values = row.asList();

        PitzerInteractionParamAttribs attribs;
        attribs.model = extractCorrectionModel(row, values);
        attribs.formulas = extractChemicalFormulasBinary(row, values);
        attribs.coefficients = extractCoefficientsBinary(row, values);
        return attribs;
    };

    auto decodeInteractionParamAttribsTernary = [&](Data const& row) -> PitzerInteractionParamAttribs
    {
        errorifnot(row.isList(), "While decoding the data shown next, I was expecting a YAML/JSON list but got instead:\n", row.dumpYaml());
        const auto values = row.asList();

        PitzerInteractionParamAttribs attribs;
        attribs.model = extractCorrectionModel(row, values);
        attribs.formulas = extractChemicalFormulasTernary(row, values);
        attribs.coefficients = extractCoefficientsTernary(row, values);
        return attribs;
    };

    auto decodeAlphaParamAttribs = [](Data const& row) -> PitzerAlphaParamAttribs
    {
        errorifnot(row.isList(), "While decoding the data shown next, I was expecting a YAML/JSON list but got instead:\n", row.dumpYaml());
        const auto values = row.asList();

        errorifnot(values.size() == 3, "Expecting a list with three values such as [Na+, Cl-, 2.3] for specification of Pitzer's alpha parameter.");
        errorifnot(values[0].isString(), "Expecting the 1st entry in the following list to be a valid chemical formula:\n", row.dumpYaml());
        errorifnot(values[1].isString(), "Expecting the 2nd entry in the following list to be a valid chemical formula:\n", row.dumpYaml());
        errorifnot(values[2].isParam(), "Expecting the 3rd entry in the following list to be a valid number:\n", row.dumpYaml());

        PitzerAlphaParamAttribs attribs;
        attribs.formula1 = values[0].asString();
        attribs.formula2 = values[1].asString();
        attribs.value = values[2].asParam();

        return attribs;
    };

    if(data.exists("beta0")) {
        errorifnot(data["beta0"].isList(), "Expecting Pitzer data for parameter beta0 to be provided as a YAML/JSON list but got this instead:\n", data["beta0"].dumpYaml());
        for(auto const& row : data["beta0"].asList())
            obj.beta0.push_back( decodeInteractionParamAttribsBinary(row) );
    }

    if(data.exists("beta1")) {
        errorifnot(data["beta1"].isList(), "Expecting Pitzer data for parameter beta1 to be provided as a YAML/JSON list but got this instead:\n", data["beta1"].dumpYaml());
        for(auto const& row : data["beta1"].asList())
            obj.beta1.push_back( decodeInteractionParamAttribsBinary(row) );
    }

    if(data.exists("beta2")) {
        errorifnot(data["beta2"].isList(), "Expecting Pitzer data for parameter beta2 to be provided as a YAML/JSON list but got this instead:\n", data["beta2"].dumpYaml());
        for(auto const& row : data["beta2"].asList())
            obj.beta2.push_back( decodeInteractionParamAttribsBinary(row) );
    }

    if(data.exists("Cphi")) {
        errorifnot(data["Cphi"].isList(), "Expecting Pitzer data for parameter Cphi to be provided as a YAML/JSON list but got this instead:\n", data["Cphi"].dumpYaml());
        for(auto const& row : data["Cphi"].asList())
            obj.Cphi.push_back( decodeInteractionParamAttribsBinary(row) );
    }

    if(data.exists("theta")) {
        errorifnot(data["theta"].isList(), "Expecting Pitzer data for parameter theta to be provided as a YAML/JSON list but got this instead:\n", data["theta"].dumpYaml());
        for(auto const& row : data["theta"].asList())
            obj.theta.push_back( decodeInteractionParamAttribsBinary(row) );
    }

    if(data.exists("psi")) {
        errorifnot(data["psi"].isList(), "Expecting Pitzer data for parameter psi to be provided as a YAML/JSON list but got this instead:\n", data["psi"].dumpYaml());
        for(auto const& row : data["psi"].asList())
            obj.psi.push_back( decodeInteractionParamAttribsTernary(row) );
    }

    if(data.exists("lambda")) {
        errorifnot(data["lambda"].isList(), "Expecting Pitzer data for parameter lambda to be provided as a YAML/JSON list but got this instead:\n", data["lambda"].dumpYaml());
        for(auto const& row : data["lambda"].asList())
            obj.lambda.push_back( decodeInteractionParamAttribsBinary(row) );
    }

    if(data.exists("zeta")) {
        errorifnot(data["zeta"].isList(), "Expecting Pitzer data for parameter zeta to be provided as a YAML/JSON list but got this instead:\n", data["zeta"].dumpYaml());
        for(auto const& row : data["zeta"].asList())
            obj.zeta.push_back( decodeInteractionParamAttribsTernary(row) );
    }

    if(data.exists("mu")) {
        errorifnot(data["mu"].isList(), "Expecting Pitzer data for parameter mu to be provided as a YAML/JSON list but got this instead:\n", data["mu"].dumpYaml());
        for(auto const& row : data["mu"].asList())
            obj.mu.push_back( decodeInteractionParamAttribsTernary(row) );
    }

    if(data.exists("eta")) {
        errorifnot(data["eta"].isList(), "Expecting Pitzer data for parameter eta to be provided as a YAML/JSON list but got this instead:\n", data["eta"].dumpYaml());
        for(auto const& row : data["eta"].asList())
            obj.eta.push_back( decodeInteractionParamAttribsTernary(row) );
    }

    if(data.exists("alpha1")) {
        errorifnot(data["alpha1"].isList(), "Expecting Pitzer data for parameter alpha1 to be provided as a YAML/JSON list but got this instead:\n", data["alpha1"].dumpYaml());
        for(auto const& row : data["alpha1"].asList())
            obj.alpha1.push_back( decodeAlphaParamAttribs(row) );
    }

    if(data.exists("alpha2")) {
        errorifnot(data["alpha2"].isList(), "Expecting Pitzer data for parameter alpha2 to be provided as a YAML/JSON list but got this instead:\n", data["alpha2"].dumpYaml());
        for(auto const& row : data["alpha2"].asList())
            obj.alpha2.push_back( decodeAlphaParamAttribs(row) );
    }
}

} // namespace Reaktoro
