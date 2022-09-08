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

#include "Models.YAML.hpp"

// Reaktoro includes
#include <Reaktoro/Models/ReactionRateModels/ReactionRateModelPalandriKharaka.hpp>
#include <Reaktoro/Models/ReactionThermoModels/ReactionThermoModelConstLgK.hpp>
#include <Reaktoro/Models/ReactionThermoModels/ReactionThermoModelGemsLgK.hpp>
#include <Reaktoro/Models/ReactionThermoModels/ReactionThermoModelPhreeqcLgK.hpp>
#include <Reaktoro/Models/ReactionThermoModels/ReactionThermoModelVantHoff.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelConstant.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelHKF.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelHollandPowell.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelInterpolation.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelMaierKelley.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelMineralHKF.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelNasa.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelWaterHKF.hpp>
#include <Reaktoro/Models/StandardVolumeModels/StandardVolumeModelConstant.hpp>
#include <Reaktoro/Serialization/Common.YAML.hpp>
#include <Reaktoro/Serialization/Core.YAML.hpp>

namespace Reaktoro {

//======================================================================
// ReactionThermoModelParams Types
//======================================================================

REAKTORO_YAML_ENCODE_DEFINE(ReactionThermoModelParamsConstLgK)
{
    node["lgKr"] = obj.lgKr;
    node["Pr"]   = obj.Pr;
}

REAKTORO_YAML_DECODE_DEFINE(ReactionThermoModelParamsConstLgK)
{
    static const ReactionThermoModelParamsConstLgK defaultstate;
    node.at("lgKr").to(obj.lgKr);
    node["Pr"].to(obj.Pr, defaultstate.Pr);
}

//----------------------------------------------------------------------

REAKTORO_YAML_ENCODE_DEFINE(ReactionThermoModelParamsGemsLgK)
{
    node["A0"] = obj.A0;
    node["A1"] = obj.A1;
    node["A2"] = obj.A2;
    node["A3"] = obj.A3;
    node["A4"] = obj.A4;
    node["A5"] = obj.A5;
    node["A6"] = obj.A6;
    node["Pr"] = obj.Pr;
}

REAKTORO_YAML_DECODE_DEFINE(ReactionThermoModelParamsGemsLgK)
{
    static const ReactionThermoModelParamsGemsLgK defaultstate;
    node.at("A0").to(obj.A0);
    node.at("A1").to(obj.A1);
    node.at("A2").to(obj.A2);
    node.at("A3").to(obj.A3);
    node.at("A4").to(obj.A4);
    node.at("A5").to(obj.A5);
    node.at("A6").to(obj.A6);
    node["Pr"].to(obj.Pr, defaultstate.Pr);
}

//----------------------------------------------------------------------

REAKTORO_YAML_ENCODE_DEFINE(ReactionThermoModelParamsPhreeqcLgK)
{
    node["A1"] = obj.A1;
    node["A2"] = obj.A2;
    node["A3"] = obj.A3;
    node["A4"] = obj.A4;
    node["A5"] = obj.A5;
    node["A6"] = obj.A6;
    node["Pr"] = obj.Pr;
}

REAKTORO_YAML_DECODE_DEFINE(ReactionThermoModelParamsPhreeqcLgK)
{
    static const ReactionThermoModelParamsPhreeqcLgK defaultstate;
    node.at("A1").to(obj.A1);
    node.at("A2").to(obj.A2);
    node.at("A3").to(obj.A3);
    node.at("A4").to(obj.A4);
    node.at("A5").to(obj.A5);
    node.at("A6").to(obj.A6);
    node["Pr"].to(obj.Pr, defaultstate.Pr);
}

//----------------------------------------------------------------------

REAKTORO_YAML_ENCODE_DEFINE(ReactionThermoModelParamsVantHoff)
{
    node["lgKr"] = obj.lgKr;
    node["dHr"]  = obj.dHr;
    node["Tr"]   = obj.Tr;
    node["Pr"]   = obj.Pr;
}

REAKTORO_YAML_DECODE_DEFINE(ReactionThermoModelParamsVantHoff)
{
    static const ReactionThermoModelParamsVantHoff defaultstate;
    node.at("lgKr").to(obj.lgKr);
    node.at("dHr").to(obj.dHr);
    node["Tr"].to(obj.Tr, defaultstate.Tr);
    node["Pr"].to(obj.Pr, defaultstate.Pr);
}

//======================================================================
// StandardThermoModelParams Types
//======================================================================

REAKTORO_YAML_ENCODE_DEFINE(StandardThermoModelParamsConstant)
{
    node["G0"]  = obj.G0;
    node["H0"]  = obj.H0;
    node["V0"]  = obj.V0;
    node["VT0"] = obj.VT0;
    node["VP0"] = obj.VP0;
    node["Cp0"] = obj.Cp0;
}

REAKTORO_YAML_DECODE_DEFINE(StandardThermoModelParamsConstant)
{
    node.at("G0").to(obj.G0);
    node["H0"].to(obj.H0, 0.0);
    node["V0"].to(obj.V0, 0.0);
    node["VT0"].to(obj.VT0, 0.0);
    node["VP0"].to(obj.VP0, 0.0);
    node["Cp0"].to(obj.Cp0, 0.0);
}

//----------------------------------------------------------------------

REAKTORO_YAML_ENCODE_DEFINE(StandardThermoModelParamsHKF)
{
    node["Gf"]     = obj.Gf;
    node["Hf"]     = obj.Hf;
    node["Sr"]     = obj.Sr;
    node["a1"]     = obj.a1;
    node["a2"]     = obj.a2;
    node["a3"]     = obj.a3;
    node["a4"]     = obj.a4;
    node["c1"]     = obj.c1;
    node["c2"]     = obj.c2;
    node["wref"]   = obj.wref;
    node["charge"] = obj.charge;
    node["Tmax"]   = obj.Tmax;
}

REAKTORO_YAML_DECODE_DEFINE(StandardThermoModelParamsHKF)
{
    node.at("Gf").to(obj.Gf);
    node.at("Hf").to(obj.Hf);
    node.at("Sr").to(obj.Sr);
    node.at("a1").to(obj.a1);
    node.at("a2").to(obj.a2);
    node.at("a3").to(obj.a3);
    node.at("a4").to(obj.a4);
    node.at("c1").to(obj.c1);
    node.at("c2").to(obj.c2);
    node.at("wref").to(obj.wref);
    node.at("charge").to(obj.charge);
}

//----------------------------------------------------------------------

REAKTORO_YAML_ENCODE_DEFINE(StandardThermoModelParamsHollandPowell)
{
    node["Gf"] = obj.Gf;
    node["Hf"] = obj.Hf;
    node["Sr"] = obj.Sr;
    node["Vr"] = obj.Vr;
    node["a"]  = obj.a;
    node["b"]  = obj.b;
    node["c"]  = obj.c;
    node["d"]  = obj.d;
    if(obj.kappa0 != 0.0)
    {
        node["alpha0"]   = obj.alpha0;
        node["kappa0"]   = obj.kappa0;
        node["kappa0p"]  = obj.kappa0p;
        node["kappa0pp"] = obj.kappa0pp;
        node["numatoms"] = obj.numatoms;
    }
    node["Tmax"] = obj.Tmax;
}

REAKTORO_YAML_DECODE_DEFINE(StandardThermoModelParamsHollandPowell)
{
    node.at("Gf").to(obj.Gf);
    node.at("Hf").to(obj.Hf);
    node.at("Sr").to(obj.Sr);
    node.at("Vr").to(obj.Vr);
    node.at("a").to(obj.a);
    node.at("b").to(obj.b);
    node.at("c").to(obj.c);
    node.at("d").to(obj.d);
    if(node["kappa0"].IsDefined())
    {
        node.at("alpha0").to(obj.alpha0);
        node.at("kappa0").to(obj.kappa0);
        node.at("kappa0p").to(obj.kappa0p);
        node.at("kappa0pp").to(obj.kappa0pp);
        node.at("numatoms").to(obj.numatoms);
    }
    node.at("Tmax").to(obj.Tmax);
}

//----------------------------------------------------------------------

REAKTORO_YAML_ENCODE_DEFINE(StandardThermoModelParamsInterpolation)
{
    node["Temperatures"] = obj.temperatures;
    node["Pressures"]    = obj.pressures;
    node["G0"]           = obj.G0;
    node["H0"]           = obj.H0;
    node["V0"]           = obj.V0;
    node["VT0"]          = obj.VT0;
    node["VP0"]          = obj.VP0;
    node["Cp0"]          = obj.Cp0;
}

REAKTORO_YAML_DECODE_DEFINE(StandardThermoModelParamsInterpolation)
{
    node.at("Temperatures").to(obj.temperatures);
    node.at("Pressures").to(obj.pressures);
    if(node["G0"]) obj.G0 = node["G0"].as<Vec<Vec<double>>>();
    if(node["H0"]) obj.H0 = node["H0"].as<Vec<Vec<double>>>();
    if(node["V0"]) obj.V0 = node["V0"].as<Vec<Vec<double>>>();
    if(node["VT0"]) obj.VT0 = node["VT0"].as<Vec<Vec<double>>>();
    if(node["VP0"]) obj.VP0 = node["VP0"].as<Vec<Vec<double>>>();
    if(node["Cp0"]) obj.Cp0 = node["Cp0"].as<Vec<Vec<double>>>();
}

//----------------------------------------------------------------------

REAKTORO_YAML_ENCODE_DEFINE(StandardThermoModelParamsMaierKelley)
{
    node["Gf"]   = obj.Gf;
    node["Hf"]   = obj.Hf;
    node["Sr"]   = obj.Sr;
    node["Vr"]   = obj.Vr;
    node["a"]    = obj.a;
    node["b"]    = obj.b;
    node["c"]    = obj.c;
    node["Tmax"] = obj.Tmax;
}

REAKTORO_YAML_DECODE_DEFINE(StandardThermoModelParamsMaierKelley)
{
    node.at("Gf").to(obj.Gf);
    node.at("Hf").to(obj.Hf);
    node.at("Sr").to(obj.Sr);
    node.at("Vr").to(obj.Vr);
    node.at("a").to(obj.a);
    node.at("b").to(obj.b);
    node.at("c").to(obj.c);
    node.at("Tmax").to(obj.Tmax);
}

//----------------------------------------------------------------------

REAKTORO_YAML_ENCODE_DEFINE(StandardThermoModelParamsMineralHKF)
{
    node["Gf"]     = obj.Gf;
    node["Hf"]     = obj.Hf;
    node["Sr"]     = obj.Sr;
    node["Vr"]     = obj.Vr;
    node["ntr"]    = obj.ntr;
    node["a"]      = obj.a;
    node["b"]      = obj.b;
    node["c"]      = obj.c;
    node["Ttr"]    = obj.Ttr;
    node["Htr"]    = obj.Htr;
    node["Vtr"]    = obj.Vtr;
    node["dPdTtr"] = obj.dPdTtr;
    node["Tmax"]   = obj.Tmax;
}

REAKTORO_YAML_DECODE_DEFINE(StandardThermoModelParamsMineralHKF)
{
    node.at("Gf").to(obj.Gf);
    node.at("Hf").to(obj.Hf);
    node.at("Sr").to(obj.Sr);
    node.at("Vr").to(obj.Vr);
    node.at("ntr").to(obj.ntr);
    node.at("a").to(obj.a);
    node.at("b").to(obj.b);
    node.at("c").to(obj.c);
    node.at("Ttr").to(obj.Ttr);
    node.at("Htr").to(obj.Htr);
    node.at("Vtr").to(obj.Vtr);
    node.at("dPdTtr").to(obj.dPdTtr);
    node.at("Tmax").to(obj.Tmax);
}

//----------------------------------------------------------------------

REAKTORO_YAML_ENCODE_DECLARE(StandardThermoModelParamsNasa::Polynomial);
REAKTORO_YAML_DECODE_DECLARE(StandardThermoModelParamsNasa::Polynomial);

REAKTORO_YAML_ENCODE_DEFINE(StandardThermoModelParamsNasa::Polynomial)
{
    node["Tmin"]  = obj.Tmin;
    node["Tmax"]  = obj.Tmax;
    node["Label"] = obj.label;
    node["State"] = obj.state;
    node["a1"]    = obj.a1;
    node["a2"]    = obj.a2;
    node["a3"]    = obj.a3;
    node["a4"]    = obj.a4;
    node["a5"]    = obj.a5;
    node["a6"]    = obj.a6;
    node["a7"]    = obj.a7;
    node["b1"]    = obj.b1;
    node["b2"]    = obj.b2;
}

REAKTORO_YAML_DECODE_DEFINE(StandardThermoModelParamsNasa::Polynomial)
{
    node.at("Tmin").to(obj.Tmin);
    node.at("Tmax").to(obj.Tmax);
    node.at("Label").to(obj.label);
    node.at("State").to(obj.state);
    node.at("a1").to(obj.a1);
    node.at("a2").to(obj.a2);
    node.at("a3").to(obj.a3);
    node.at("a4").to(obj.a4);
    node.at("a5").to(obj.a5);
    node.at("a6").to(obj.a6);
    node.at("a7").to(obj.a7);
    node.at("b1").to(obj.b1);
    node.at("b2").to(obj.b2);
}

//----------------------------------------------------------------------

REAKTORO_YAML_ENCODE_DEFINE(StandardThermoModelParamsNasa)
{
    if(obj.polynomials.size())
    {
        node["dHf"] = obj.dHf;
        node["dH0"] = obj.dH0;
        node["Polynomials"] = obj.polynomials;
    }
    else
    {
        node["H0"] = obj.H0;
        node["T0"] = obj.T0;
    }
}

REAKTORO_YAML_DECODE_DEFINE(StandardThermoModelParamsNasa)
{
    if(node["Polynomials"])
    {
        node.at("dHf").to(obj.dHf);
        node.at("dH0").to(obj.dH0);
        node.at("Polynomials").to(obj.polynomials);
    }
    else
    {
        node.at("H0").to(obj.H0);
        node.at("T0").to(obj.T0);
    }
}

//----------------------------------------------------------------------

REAKTORO_YAML_ENCODE_DEFINE(StandardThermoModelParamsWaterHKF)
{
    node["Ttr"] = obj.Ttr;
    node["Str"] = obj.Str;
    node["Gtr"] = obj.Gtr;
    node["Htr"] = obj.Htr;
}

REAKTORO_YAML_DECODE_DEFINE(StandardThermoModelParamsWaterHKF)
{
    node.at("Ttr").to(obj.Ttr);
    node.at("Str").to(obj.Str);
    node.at("Gtr").to(obj.Gtr);
    node.at("Htr").to(obj.Htr);
}

//======================================================================
// StandardVolumeModelParams Types
//======================================================================

REAKTORO_YAML_ENCODE_DEFINE(StandardVolumeModelParamsConstant)
{
    node["V0"] = obj.V0;
}

REAKTORO_YAML_DECODE_DEFINE(StandardVolumeModelParamsConstant)
{
    node.at("V0").to(obj.V0);
}

//======================================================================
// ReactionRateModelParams Types
//======================================================================

REAKTORO_YAML_ENCODE_DECLARE(ReactionRateModelParamsPalandriKharaka::Mechanism);
REAKTORO_YAML_DECODE_DECLARE(ReactionRateModelParamsPalandriKharaka::Mechanism);

REAKTORO_YAML_ENCODE_DEFINE(ReactionRateModelParamsPalandriKharaka::Mechanism)
{
    node["lgk"] = obj.lgk;
    node["E"] = obj.E;
    node.appendIfNotDefault("p", obj.p, 1.0);
    node.appendIfNotDefault("q", obj.q, 1.0);
    for(auto catalyst : obj.catalysts)
    {
        const auto key = catalyst.property + "(" + catalyst.formula + ")";
        node[key] = catalyst.power;
    }
}

REAKTORO_YAML_DECODE_DEFINE(ReactionRateModelParamsPalandriKharaka::Mechanism)
{
    node.at("lgk").to(obj.lgk);
    node.at("E").to(obj.E);
    node.copyOptionalChildValueTo("p", obj.p, Param{1.0});
    node.copyOptionalChildValueTo("q", obj.q, Param{1.0});

    // Collect all catalyst properties and their power
    for(auto child : node)
    {
        const auto key = child.first.as<String>();
        if(!oneof(key[0], 'a', 'P'))
            continue;
        errorif(key.size() <= 3, "Expecting a chemical formula inside `a()` or `P()`, such as `a(H+)`, `P(CO2)`.");
        errorif(key[1] != '(' || key.back() != ')', "Expecting ( and ) as in `a(H+)`, `a(Fe+3)`, `P(CO2)`.");
        const auto formula = key.substr(2, key.size() - 3); // exclude first two chars and last
        const auto property = key.substr(0, 1);
        const auto power = child.second.as<double>();
        obj.catalysts.push_back({ formula, property, power });
    }
}

REAKTORO_YAML_ENCODE_DEFINE(ReactionRateModelParamsPalandriKharaka)
{
    node["Mineral"] = join(obj.names, " ");
    for(auto mechanism : obj.mechanisms)
        node["Mechanisms"][mechanism.name] = mechanism;
}

REAKTORO_YAML_DECODE_DEFINE(ReactionRateModelParamsPalandriKharaka)
{
    obj.names = split(node.at("Mineral"), " ");
    for(auto mechanism : node["Mechanisms"])
    {
        obj.mechanisms.push_back(mechanism.second.as<ReactionRateModelParamsPalandriKharaka::Mechanism>());
        obj.mechanisms.back().name = mechanism.first.as<String>();
    }
}

//----------------------------------------------------------------------

} // namespace YAML
