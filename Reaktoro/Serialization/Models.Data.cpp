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

#include "Models.Data.hpp"

// Reaktoro includes
#include <Reaktoro/Models/ReactionRateModels/ReactionRateModelPalandriKharaka.hpp>
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelConstLgK.hpp>
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelGemsLgK.hpp>
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelPhreeqcLgK.hpp>
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelVantHoff.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelConstant.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelHKF.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelHollandPowell.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelInterpolation.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelMaierKelley.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelMineralHKF.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelNasa.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelWaterHKF.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardVolumeModelConstant.hpp>

namespace Reaktoro {

//======================================================================
// StandardThermoModelParams Types
//======================================================================

REAKTORO_DATA_ENCODE_DEFINE(StandardThermoModelParamsConstant)
{
    data.at("G0") = obj.G0;
    data.at("H0") = obj.H0;
    data.at("V0") = obj.V0;
    data.at("VT0") = obj.VT0;
    data.at("VP0") = obj.VP0;
    data.at("Cp0") = obj.Cp0;
}

REAKTORO_DATA_DECODE_DEFINE(StandardThermoModelParamsConstant)
{
    data.at("G0").to(obj.G0);
    data.at("H0").to(obj.H0);
    data.at("V0").to(obj.V0);
    data.at("VT0").to(obj.VT0);
    data.at("VP0").to(obj.VP0);
    data.at("Cp0").to(obj.Cp0);
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DEFINE(StandardThermoModelParamsHKF)
{
    data.at("Gf") = obj.Gf;
    data.at("Hf") = obj.Hf;
    data.at("Sr") = obj.Sr;
    data.at("a1") = obj.a1;
    data.at("a2") = obj.a2;
    data.at("a3") = obj.a3;
    data.at("a4") = obj.a4;
    data.at("c1") = obj.c1;
    data.at("c2") = obj.c2;
    data.at("wref") = obj.wref;
    data.at("charge") = obj.charge;
    data.at("Tmax") = obj.Tmax;
}

REAKTORO_DATA_DECODE_DEFINE(StandardThermoModelParamsHKF)
{
    data["Gf"].to(obj.Gf);
    data["Hf"].to(obj.Hf);
    data["Sr"].to(obj.Sr);
    data["a1"].to(obj.a1);
    data["a2"].to(obj.a2);
    data["a3"].to(obj.a3);
    data["a4"].to(obj.a4);
    data["c1"].to(obj.c1);
    data["c2"].to(obj.c2);
    data["wref"].to(obj.wref);
    data["charge"].to(obj.charge);
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DEFINE(StandardThermoModelParamsHollandPowell)
{
    data.at("Gf") = obj.Gf;
    data.at("Hf") = obj.Hf;
    data.at("Sr") = obj.Sr;
    data.at("Vr") = obj.Vr;
    data.at("a")  = obj.a;
    data.at("b")  = obj.b;
    data.at("c")  = obj.c;
    data.at("d")  = obj.d;
    if(obj.kappa0 != 0.0)
    {
        data.at("alpha0")   = obj.alpha0;
        data.at("kappa0")   = obj.kappa0;
        data.at("kappa0p")  = obj.kappa0p;
        data.at("kappa0pp") = obj.kappa0pp;
        data.at("numatoms") = obj.numatoms;
    }
    data.at("Tmax") = obj.Tmax;
}

REAKTORO_DATA_DECODE_DEFINE(StandardThermoModelParamsHollandPowell)
{
    data.at("Gf").to(obj.Gf);
    data.at("Hf").to(obj.Hf);
    data.at("Sr").to(obj.Sr);
    data.at("Vr").to(obj.Vr);
    data.at("a").to(obj.a);
    data.at("b").to(obj.b);
    data.at("c").to(obj.c);
    data.at("d").to(obj.d);
    if(data.exists("kappa0"))
    {
        data.at("alpha0").to(obj.alpha0);
        data.at("kappa0").to(obj.kappa0);
        data.at("kappa0p").to(obj.kappa0p);
        data.at("kappa0pp").to(obj.kappa0pp);
        data.at("numatoms").to(obj.numatoms);
    }
    data.at("Tmax").to(obj.Tmax);
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DEFINE(StandardThermoModelParamsInterpolation)
{
    data.at("Temperatures") = obj.temperatures;
    data.at("Pressures")    = obj.pressures;
    data.at("G0")           = obj.G0;
    data.at("H0")           = obj.H0;
    data.at("V0")           = obj.V0;
    data.at("VT0")          = obj.VT0;
    data.at("VP0")          = obj.VP0;
    data.at("Cp0")          = obj.Cp0;
}

REAKTORO_DATA_DECODE_DEFINE(StandardThermoModelParamsInterpolation)
{
    data.at("Temperatures").to(obj.temperatures);
    data.at("Pressures").to(obj.pressures);
    data.at("G0").to(obj.G0);
    if(data.exists("H0")) data.at("H0").to(obj.H0);
    if(data.exists("V0")) data.at("V0").to(obj.V0);
    if(data.exists("VT0")) data.at("VT0").to(obj.VT0);
    if(data.exists("VP0")) data.at("VP0").to(obj.VP0);
    if(data.exists("Cp0")) data.at("Cp0").to(obj.Cp0);
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DEFINE(StandardThermoModelParamsMaierKelley)
{
    data.at("Gf")   = obj.Gf;
    data.at("Hf")   = obj.Hf;
    data.at("Sr")   = obj.Sr;
    data.at("Vr")   = obj.Vr;
    data.at("a")    = obj.a;
    data.at("b")    = obj.b;
    data.at("c")    = obj.c;
    data.at("Tmax") = obj.Tmax;
}

REAKTORO_DATA_DECODE_DEFINE(StandardThermoModelParamsMaierKelley)
{
    data.at("Gf").to(obj.Gf);
    data.at("Hf").to(obj.Hf);
    data.at("Sr").to(obj.Sr);
    data.at("Vr").to(obj.Vr);
    data.at("a").to(obj.a);
    data.at("b").to(obj.b);
    data.at("c").to(obj.c);
    data.at("Tmax").to(obj.Tmax);
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DEFINE(StandardThermoModelParamsMineralHKF)
{
    data.at("Gf")     = obj.Gf;
    data.at("Hf")     = obj.Hf;
    data.at("Sr")     = obj.Sr;
    data.at("Vr")     = obj.Vr;
    data.at("ntr")    = obj.ntr;
    data.at("a")      = obj.a;
    data.at("b")      = obj.b;
    data.at("c")      = obj.c;
    data.at("Ttr")    = obj.Ttr;
    data.at("Htr")    = obj.Htr;
    data.at("Vtr")    = obj.Vtr;
    data.at("dPdTtr") = obj.dPdTtr;
    data.at("Tmax")   = obj.Tmax;
}

REAKTORO_DATA_DECODE_DEFINE(StandardThermoModelParamsMineralHKF)
{
    data.at("Gf").to(obj.Gf);
    data.at("Hf").to(obj.Hf);
    data.at("Sr").to(obj.Sr);
    data.at("Vr").to(obj.Vr);
    data.at("ntr").to(obj.ntr);
    data.at("a").to(obj.a);
    data.at("b").to(obj.b);
    data.at("c").to(obj.c);
    data.at("Ttr").to(obj.Ttr);
    data.at("Htr").to(obj.Htr);
    data.at("Vtr").to(obj.Vtr);
    data.at("dPdTtr").to(obj.dPdTtr);
    data.at("Tmax").to(obj.Tmax);
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DECLARE(StandardThermoModelParamsNasa::Polynomial);
REAKTORO_DATA_DECODE_DECLARE(StandardThermoModelParamsNasa::Polynomial);

REAKTORO_DATA_ENCODE_DEFINE(StandardThermoModelParamsNasa::Polynomial)
{
    data.at("Tmin")  = obj.Tmin;
    data.at("Tmax")  = obj.Tmax;
    data.at("Label") = obj.label;
    data.at("State") = obj.state;
    data.at("a1")    = obj.a1;
    data.at("a2")    = obj.a2;
    data.at("a3")    = obj.a3;
    data.at("a4")    = obj.a4;
    data.at("a5")    = obj.a5;
    data.at("a6")    = obj.a6;
    data.at("a7")    = obj.a7;
    data.at("b1")    = obj.b1;
    data.at("b2")    = obj.b2;
}

REAKTORO_DATA_DECODE_DEFINE(StandardThermoModelParamsNasa::Polynomial)
{
    data.at("Tmin").to(obj.Tmin);
    data.at("Tmax").to(obj.Tmax);
    data.at("Label").to(obj.label);
    data.at("State").to(obj.state);
    data.at("a1").to(obj.a1);
    data.at("a2").to(obj.a2);
    data.at("a3").to(obj.a3);
    data.at("a4").to(obj.a4);
    data.at("a5").to(obj.a5);
    data.at("a6").to(obj.a6);
    data.at("a7").to(obj.a7);
    data.at("b1").to(obj.b1);
    data.at("b2").to(obj.b2);
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DEFINE(StandardThermoModelParamsNasa)
{
    if(obj.polynomials.size())
    {
        data.at("dHf") = obj.dHf;
        data.at("dH0") = obj.dH0;
        data.at("Polynomials") = obj.polynomials;
    }
    else
    {
        data.at("H0") = obj.H0;
        data.at("T0") = obj.T0;
    }
}

REAKTORO_DATA_DECODE_DEFINE(StandardThermoModelParamsNasa)
{
    if(data.exists("Polynomials"))
    {
        data.at("dHf").to(obj.dHf);
        data.at("dH0").to(obj.dH0);
        data.at("Polynomials").to(obj.polynomials);
    }
    else
    {
        data.at("H0").to(obj.H0);
        data.at("T0").to(obj.T0);
    }
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DEFINE(StandardThermoModelParamsWaterHKF)
{
    data.at("Ttr") = obj.Ttr;
    data.at("Str") = obj.Str;
    data.at("Gtr") = obj.Gtr;
    data.at("Htr") = obj.Htr;
}

REAKTORO_DATA_DECODE_DEFINE(StandardThermoModelParamsWaterHKF)
{
    data.at("Ttr").to(obj.Ttr);
    data.at("Str").to(obj.Str);
    data.at("Gtr").to(obj.Gtr);
    data.at("Htr").to(obj.Htr);
}

//======================================================================
// ReactionStandardThermoModelParams Types
//======================================================================

REAKTORO_DATA_ENCODE_DEFINE(ReactionStandardThermoModelParamsConstLgK)
{
    data.at("lgKr") = obj.lgKr;
    data.at("Pr")   = obj.Pr;
}

REAKTORO_DATA_DECODE_DEFINE(ReactionStandardThermoModelParamsConstLgK)
{
    data.at("lgKr").to(obj.lgKr);
    if(data.exists("Pr")) data.at("Pr").to(obj.Pr);
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DEFINE(ReactionStandardThermoModelParamsGemsLgK)
{
    data.at("A0") = obj.A0;
    data.at("A1") = obj.A1;
    data.at("A2") = obj.A2;
    data.at("A3") = obj.A3;
    data.at("A4") = obj.A4;
    data.at("A5") = obj.A5;
    data.at("A6") = obj.A6;
    data.at("Pr") = obj.Pr;
}

REAKTORO_DATA_DECODE_DEFINE(ReactionStandardThermoModelParamsGemsLgK)
{
    data.at("A0").to(obj.A0);
    data.at("A1").to(obj.A1);
    data.at("A2").to(obj.A2);
    data.at("A3").to(obj.A3);
    data.at("A4").to(obj.A4);
    data.at("A5").to(obj.A5);
    data.at("A6").to(obj.A6);
    if(data.exists("Pr")) data.at("Pr").to(obj.Pr);
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DEFINE(ReactionStandardThermoModelParamsPhreeqcLgK)
{
    data.at("A1") = obj.A1;
    data.at("A2") = obj.A2;
    data.at("A3") = obj.A3;
    data.at("A4") = obj.A4;
    data.at("A5") = obj.A5;
    data.at("A6") = obj.A6;
    data.at("Pr") = obj.Pr;
}

REAKTORO_DATA_DECODE_DEFINE(ReactionStandardThermoModelParamsPhreeqcLgK)
{
    data.at("A1").to(obj.A1);
    data.at("A2").to(obj.A2);
    data.at("A3").to(obj.A3);
    data.at("A4").to(obj.A4);
    data.at("A5").to(obj.A5);
    data.at("A6").to(obj.A6);
    data.at("Pr").to(obj.Pr);
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DEFINE(ReactionStandardThermoModelParamsVantHoff)
{
    data.at("lgKr") = obj.lgKr;
    data.at("dHr")  = obj.dHr;
    if(data.exists("Pr")) data.at("Tr") = obj.Tr;
    if(data.exists("Pr")) data.at("Pr") = obj.Pr;
}

REAKTORO_DATA_DECODE_DEFINE(ReactionStandardThermoModelParamsVantHoff)
{
    data.at("lgKr").to(obj.lgKr);
    data.at("dHr").to(obj.dHr);
    if(data.exists("Tr")) data.at("Tr").to(obj.Tr);
    if(data.exists("Pr")) data.at("Pr").to(obj.Pr);
}

//======================================================================
// StandardVolumeModelParams Types
//======================================================================

REAKTORO_DATA_ENCODE_DEFINE(StandardVolumeModelParamsConstant)
{
    data.at("V0") = obj.V0;
}

REAKTORO_DATA_DECODE_DEFINE(StandardVolumeModelParamsConstant)
{
    data.at("V0").to(obj.V0);
}

//======================================================================
// ReactionRateModelParams Types
//======================================================================

REAKTORO_DATA_ENCODE_DECLARE(ReactionRateModelParamsPalandriKharaka::Mechanism);
REAKTORO_DATA_DECODE_DECLARE(ReactionRateModelParamsPalandriKharaka::Mechanism);

REAKTORO_DATA_ENCODE_DEFINE(ReactionRateModelParamsPalandriKharaka::Mechanism)
{
    data.at("lgk") = obj.lgk;
    data.at("E") = obj.E;
    data.at("p") = obj.p;
    data.at("q") = obj.q;
    for(auto catalyst : obj.catalysts)
    {
        const auto key = catalyst.property + "(" + catalyst.formula + ")";
        data.at(key) = catalyst.power;
    }
}

REAKTORO_DATA_DECODE_DEFINE(ReactionRateModelParamsPalandriKharaka::Mechanism)
{
    data.at("lgk").to(obj.lgk);
    data.at("E").to(obj.E);
    if(data.exists("p")) data.at("p").to(obj.p);
    if(data.exists("q")) data.at("q").to(obj.q);

    // Collect all catalyst properties and their power
    for(auto const& [key, value] : data.asDict())
    {
        if(!oneof(key[0], 'a', 'P'))
            continue;
        errorif(key.size() <= 3, "Expecting a chemical formula inside `a()` or `P()`, such as `a(H+)`, `P(CO2)`.");
        errorif(key[1] != '(' || key.back() != ')', "Expecting ( and ) as in `a(H+)`, `a(Fe+3)`, `P(CO2)`.");
        const auto formula = key.substr(2, key.size() - 3); // exclude first two chars and last
        const auto property = key.substr(0, 1);
        const auto power = value.asFloat();
        obj.catalysts.push_back({ formula, property, power });
    }
}

REAKTORO_DATA_ENCODE_DEFINE(ReactionRateModelParamsPalandriKharaka)
{
    data.at("Mineral") = join(obj.names, " ");
    for(auto mechanism : obj.mechanisms)
        data.at("Mechanisms").at(mechanism.name) = mechanism;
}

REAKTORO_DATA_DECODE_DEFINE(ReactionRateModelParamsPalandriKharaka)
{
    obj.names = split(data.at("Mineral"), " ");
    for(auto const& [name, mechanism] : data.at("Mechanisms").asDict())
    {
        obj.mechanisms.push_back(mechanism.as<ReactionRateModelParamsPalandriKharaka::Mechanism>());
        obj.mechanisms.back().name = name;
    }
}

//----------------------------------------------------------------------

} // namespace YAML
