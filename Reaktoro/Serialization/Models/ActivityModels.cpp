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

REAKTORO_DATA_ENCODE_DECLARE(ActivityModelParamsPitzer::CorrectionModel);
REAKTORO_DATA_DECODE_DECLARE(ActivityModelParamsPitzer::CorrectionModel);

REAKTORO_DATA_ENCODE_DEFINE(ActivityModelParamsPitzer::CorrectionModel)
{
    auto tostring = [](ActivityModelParamsPitzer::CorrectionModel modeltype) -> String
    {
        switch(modeltype)
        {
            case ActivityModelParamsPitzer::CorrectionModel::Constant:           return "Constant";
            case ActivityModelParamsPitzer::CorrectionModel::Phreeqc:            return "Phreeqc";
            case ActivityModelParamsPitzer::CorrectionModel::HeMorse1993:        return "HeMorse1993";
            case ActivityModelParamsPitzer::CorrectionModel::Dai2013:            return "Dai2013";
            case ActivityModelParamsPitzer::CorrectionModel::Dai2014:            return "Dai2014";
            case ActivityModelParamsPitzer::CorrectionModel::ChristovMoller2004: return "ChristovMoller2004";
            case ActivityModelParamsPitzer::CorrectionModel::Holmes1987:         return "Holmes1987";
            case ActivityModelParamsPitzer::CorrectionModel::Pitzer1984:         return "Pitzer1984";
            case ActivityModelParamsPitzer::CorrectionModel::PalabanPitzer1987:  return "PalabanPitzer1987";
            case ActivityModelParamsPitzer::CorrectionModel::Polya2001:          return "Polya2001";
            case ActivityModelParamsPitzer::CorrectionModel::LiDuan2007:         return "LiDuan2007";
        }
        errorif(true, "Unable to encode Pitzer parameter (T, P) correction model `", modeltype, "` to a string.");
        return {};
    };

    data = tostring(obj);
}

REAKTORO_DATA_DECODE_DEFINE(ActivityModelParamsPitzer::CorrectionModel)
{
    auto fromstring = [](String const& modeltype) -> ActivityModelParamsPitzer::CorrectionModel
    {
        if(modeltype == "Constant")            return ActivityModelParamsPitzer::CorrectionModel::Constant;
        if(modeltype == "Phreeqc")             return ActivityModelParamsPitzer::CorrectionModel::Phreeqc;
        if(modeltype == "HeMorse1993")         return ActivityModelParamsPitzer::CorrectionModel::HeMorse1993;
        if(modeltype == "Dai2013")             return ActivityModelParamsPitzer::CorrectionModel::Dai2013;
        if(modeltype == "Dai2014")             return ActivityModelParamsPitzer::CorrectionModel::Dai2014;
        if(modeltype == "ChristovMoller2004")  return ActivityModelParamsPitzer::CorrectionModel::ChristovMoller2004;
        if(modeltype == "Holmes1987")          return ActivityModelParamsPitzer::CorrectionModel::Holmes1987;
        if(modeltype == "Pitzer1984")          return ActivityModelParamsPitzer::CorrectionModel::Pitzer1984;
        if(modeltype == "PalabanPitzer1987")   return ActivityModelParamsPitzer::CorrectionModel::PalabanPitzer1987;
        if(modeltype == "Polya2001")           return ActivityModelParamsPitzer::CorrectionModel::Polya2001;
        if(modeltype == "LiDuan2007")          return ActivityModelParamsPitzer::CorrectionModel::LiDuan2007;

        errorif(true, "Unable to decode Pitzer parameter (T, P) correction model from string `", modeltype, "`.");
        return {};
    };

    errorifnot(data.isString(), "Expecting a string representative of a ActivityModelParamsPitzer::CorrectionModel value but got this instead:\n", data.dump());

    obj = fromstring(data.asString());
}

REAKTORO_DATA_ENCODE_DECLARE(ActivityModelParamsPitzer::InteractionParamAttribs);
REAKTORO_DATA_DECODE_DECLARE(ActivityModelParamsPitzer::InteractionParamAttribs);

REAKTORO_DATA_ENCODE_DEFINE(ActivityModelParamsPitzer::InteractionParamAttribs)
{
    data["Formulas"] = obj.formulas;
    data["CorrectionModel"] = obj.model;
    data["Parameters"] = obj.parameters;
}

REAKTORO_DATA_DECODE_DEFINE(ActivityModelParamsPitzer::InteractionParamAttribs)
{
    data.required("Formulas").to(obj.formulas);
    data.optional("CorrectionModel").to(obj.model);
    data.required("Parameters").to(obj.parameters);
}

REAKTORO_DATA_ENCODE_DEFINE(ActivityModelParamsPitzer)
{
    data["Beta0"]  = obj.beta0;
    data["Beta1"]  = obj.beta1;
    data["Beta2"]  = obj.beta2;
    data["Cphi"]   = obj.Cphi;
    data["Theta"]  = obj.theta;
    data["Psi"]    = obj.psi;
    data["Lambda"] = obj.lambda;
    data["Zeta"]   = obj.zeta;
    data["Mu"]     = obj.mu;
    data["Eta"]    = obj.eta;
    data["Alpha1"] = obj.alpha1;
    data["Alpha2"] = obj.alpha2;
}

REAKTORO_DATA_DECODE_DEFINE(ActivityModelParamsPitzer)
{
    data.optional("Beta0").to(obj.beta0);
    data.optional("Beta1").to(obj.beta1);
    data.optional("Beta2").to(obj.beta2);
    data.optional("Cphi").to(obj.Cphi);
    data.optional("Theta").to(obj.theta);
    data.optional("Psi").to(obj.psi);
    data.optional("Lambda").to(obj.lambda);
    data.optional("Zeta").to(obj.zeta);
    data.optional("Mu").to(obj.mu);
    data.optional("Eta").to(obj.eta);
    data.optional("Alpha1").to(obj.alpha1);
    data.optional("Alpha2").to(obj.alpha2);
}

} // namespace Reaktoro
