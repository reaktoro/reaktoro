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

#include "StandardThermoModels.hpp"

// Reaktoro includes
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Models/ReactionRateModels.hpp>
#include <Reaktoro/Models/StandardThermoModels.hpp>
#include <Reaktoro/Serialization/Common.hpp>
#include <Reaktoro/Serialization/Core.hpp>

namespace Reaktoro {

//======================================================================
// StandardThermoModelParams Types
//======================================================================

REAKTORO_DATA_ENCODE_DEFINE(StandardThermoModelParamsConstant)
{
    data["G0"] = obj.G0;
    data["H0"] = obj.H0;
    data["V0"] = obj.V0;
    data["VT0"] = obj.VT0;
    data["VP0"] = obj.VP0;
    data["Cp0"] = obj.Cp0;
}

REAKTORO_DATA_DECODE_DEFINE(StandardThermoModelParamsConstant)
{
    data.required("G0").to(obj.G0);
    data.optional("H0").to(obj.H0);
    data.optional("V0").to(obj.V0);
    data.optional("VT0").to(obj.VT0);
    data.optional("VP0").to(obj.VP0);
    data.optional("Cp0").to(obj.Cp0);
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DEFINE(StandardThermoModelParamsHKF)
{
    data["Gf"] = obj.Gf;
    data["Hf"] = obj.Hf;
    data["Sr"] = obj.Sr;
    data["a1"] = obj.a1;
    data["a2"] = obj.a2;
    data["a3"] = obj.a3;
    data["a4"] = obj.a4;
    data["c1"] = obj.c1;
    data["c2"] = obj.c2;
    data["wref"] = obj.wref;
    data["charge"] = obj.charge;
    data["Tmax"] = obj.Tmax;
}

REAKTORO_DATA_DECODE_DEFINE(StandardThermoModelParamsHKF)
{
    data.required("Gf").to(obj.Gf);
    data.required("Hf").to(obj.Hf);
    data.required("Sr").to(obj.Sr);
    data.required("a1").to(obj.a1);
    data.required("a2").to(obj.a2);
    data.required("a3").to(obj.a3);
    data.required("a4").to(obj.a4);
    data.required("c1").to(obj.c1);
    data.required("c2").to(obj.c2);
    data.required("wref").to(obj.wref);
    data.required("charge").to(obj.charge);
    data.optional("Tmax").to(obj.Tmax);
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DEFINE(StandardThermoModelParamsHollandPowell)
{
    data["Gf"] = obj.Gf;
    data["Hf"] = obj.Hf;
    data["Sr"] = obj.Sr;
    data["Vr"] = obj.Vr;
    data["a"]  = obj.a;
    data["b"]  = obj.b;
    data["c"]  = obj.c;
    data["d"]  = obj.d;
    if(obj.kappa0 != 0.0)
    {
        data["alpha0"]   = obj.alpha0;
        data["kappa0"]   = obj.kappa0;
        data["kappa0p"]  = obj.kappa0p;
        data["kappa0pp"] = obj.kappa0pp;
        data["numatoms"] = obj.numatoms;
    }
    data["Tmax"] = obj.Tmax;
}

REAKTORO_DATA_DECODE_DEFINE(StandardThermoModelParamsHollandPowell)
{
    data.required("Gf").to(obj.Gf);
    data.required("Hf").to(obj.Hf);
    data.required("Sr").to(obj.Sr);
    data.required("Vr").to(obj.Vr);
    data.required("a").to(obj.a);
    data.required("b").to(obj.b);
    data.required("c").to(obj.c);
    data.required("d").to(obj.d);
    if(data.exists("kappa0"))
    {
        data.required("alpha0").to(obj.alpha0);
        data.required("kappa0").to(obj.kappa0);
        data.required("kappa0p").to(obj.kappa0p);
        data.required("kappa0pp").to(obj.kappa0pp);
        data.required("numatoms").to(obj.numatoms);
    }
    data.optional("Tmax").to(obj.Tmax);
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DEFINE(StandardThermoModelParamsInterpolation)
{
    data["Temperatures"] = obj.temperatures;
    data["Pressures"]    = obj.pressures;
    data["G0"]           = obj.G0;
    data["H0"]           = obj.H0;
    data["V0"]           = obj.V0;
    data["VT0"]          = obj.VT0;
    data["VP0"]          = obj.VP0;
    data["Cp0"]          = obj.Cp0;
}

REAKTORO_DATA_DECODE_DEFINE(StandardThermoModelParamsInterpolation)
{
    data.required("Temperatures").to(obj.temperatures);
    data.required("Pressures").to(obj.pressures);
    data.required("G0").to(obj.G0);
    data.optional("H0").to(obj.H0);
    data.optional("V0").to(obj.V0);
    data.optional("VT0").to(obj.VT0);
    data.optional("VP0").to(obj.VP0);
    data.optional("Cp0").to(obj.Cp0);
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DEFINE(StandardThermoModelParamsMaierKelley)
{
    data["Gf"]   = obj.Gf;
    data["Hf"]   = obj.Hf;
    data["Sr"]   = obj.Sr;
    data["Vr"]   = obj.Vr;
    data["a"]    = obj.a;
    data["b"]    = obj.b;
    data["c"]    = obj.c;
    data["Tmax"] = obj.Tmax;
}

REAKTORO_DATA_DECODE_DEFINE(StandardThermoModelParamsMaierKelley)
{
    data.required("Gf").to(obj.Gf);
    data.required("Hf").to(obj.Hf);
    data.required("Sr").to(obj.Sr);
    data.required("Vr").to(obj.Vr);
    data.required("a").to(obj.a);
    data.required("b").to(obj.b);
    data.required("c").to(obj.c);
    data.optional("Tmax").to(obj.Tmax);
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DEFINE(StandardThermoModelParamsMineralHKF)
{
    data["Gf"]     = obj.Gf;
    data["Hf"]     = obj.Hf;
    data["Sr"]     = obj.Sr;
    data["Vr"]     = obj.Vr;
    data["ntr"]    = obj.ntr;
    data["a"]      = obj.a;
    data["b"]      = obj.b;
    data["c"]      = obj.c;
    data["Ttr"]    = obj.Ttr;
    data["Htr"]    = obj.Htr;
    data["Vtr"]    = obj.Vtr;
    data["dPdTtr"] = obj.dPdTtr;
    data["Tmax"]   = obj.Tmax;
}

REAKTORO_DATA_DECODE_DEFINE(StandardThermoModelParamsMineralHKF)
{
    data.required("Gf").to(obj.Gf);
    data.required("Hf").to(obj.Hf);
    data.required("Sr").to(obj.Sr);
    data.required("Vr").to(obj.Vr);
    data.required("ntr").to(obj.ntr);
    data.required("a").to(obj.a);
    data.required("b").to(obj.b);
    data.required("c").to(obj.c);
    data.required("Ttr").to(obj.Ttr);
    data.required("Htr").to(obj.Htr);
    data.required("Vtr").to(obj.Vtr);
    data.required("dPdTtr").to(obj.dPdTtr);
    data.optional("Tmax").to(obj.Tmax);
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DECLARE(StandardThermoModelParamsNasa::Polynomial);
REAKTORO_DATA_DECODE_DECLARE(StandardThermoModelParamsNasa::Polynomial);

REAKTORO_DATA_ENCODE_DEFINE(StandardThermoModelParamsNasa::Polynomial)
{
    data["Tmin"]  = obj.Tmin;
    data["Tmax"]  = obj.Tmax;
    data["Label"] = obj.label;
    data["State"] = obj.state;
    data["a1"]    = obj.a1;
    data["a2"]    = obj.a2;
    data["a3"]    = obj.a3;
    data["a4"]    = obj.a4;
    data["a5"]    = obj.a5;
    data["a6"]    = obj.a6;
    data["a7"]    = obj.a7;
    data["b1"]    = obj.b1;
    data["b2"]    = obj.b2;
}

REAKTORO_DATA_DECODE_DEFINE(StandardThermoModelParamsNasa::Polynomial)
{
    data.required("Tmin").to(obj.Tmin);
    data.required("Tmax").to(obj.Tmax);
    data.required("Label").to(obj.label);
    data.required("State").to(obj.state);
    data.required("a1").to(obj.a1);
    data.required("a2").to(obj.a2);
    data.required("a3").to(obj.a3);
    data.required("a4").to(obj.a4);
    data.required("a5").to(obj.a5);
    data.required("a6").to(obj.a6);
    data.required("a7").to(obj.a7);
    data.required("b1").to(obj.b1);
    data.required("b2").to(obj.b2);
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DEFINE(StandardThermoModelParamsNasa)
{
    if(obj.polynomials.size())
    {
        data["dHf"] = obj.dHf;
        data["dH0"] = obj.dH0;
        data["Polynomials"] = obj.polynomials;
    }
    else
    {
        data["H0"] = obj.H0;
        data["T0"] = obj.T0;
    }
}

REAKTORO_DATA_DECODE_DEFINE(StandardThermoModelParamsNasa)
{
    if(data.exists("Polynomials"))
    {
        data.required("dHf").to(obj.dHf);
        data.required("dH0").to(obj.dH0);
        data.required("Polynomials").to(obj.polynomials);
    }
    else
    {
        data.required("H0").to(obj.H0);
        data.required("T0").to(obj.T0);
    }
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DEFINE(StandardThermoModelParamsWaterHKF)
{
    data["Ttr"] = obj.Ttr;
    data["Str"] = obj.Str;
    data["Gtr"] = obj.Gtr;
    data["Htr"] = obj.Htr;
}

REAKTORO_DATA_DECODE_DEFINE(StandardThermoModelParamsWaterHKF)
{
    data.required("Ttr").to(obj.Ttr);
    data.required("Str").to(obj.Str);
    data.required("Gtr").to(obj.Gtr);
    data.required("Htr").to(obj.Htr);
}

//======================================================================
// ReactionStandardThermoModelParams Types
//======================================================================

REAKTORO_DATA_ENCODE_DEFINE(ReactionStandardThermoModelParamsConstLgK)
{
    data["lgKr"] = obj.lgKr;
    data["Pr"]   = obj.Pr;
}

REAKTORO_DATA_DECODE_DEFINE(ReactionStandardThermoModelParamsConstLgK)
{
    data.required("lgKr").to(obj.lgKr);
    data.optional("Pr").to(obj.Pr);
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DEFINE(ReactionStandardThermoModelParamsGemsLgK)
{
    data["A0"] = obj.A0;
    data["A1"] = obj.A1;
    data["A2"] = obj.A2;
    data["A3"] = obj.A3;
    data["A4"] = obj.A4;
    data["A5"] = obj.A5;
    data["A6"] = obj.A6;
    data["Pr"] = obj.Pr;
}

REAKTORO_DATA_DECODE_DEFINE(ReactionStandardThermoModelParamsGemsLgK)
{
    data.required("A0").to(obj.A0);
    data.required("A1").to(obj.A1);
    data.required("A2").to(obj.A2);
    data.required("A3").to(obj.A3);
    data.required("A4").to(obj.A4);
    data.required("A5").to(obj.A5);
    data.required("A6").to(obj.A6);
    data.optional("Pr").to(obj.Pr);
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DEFINE(ReactionStandardThermoModelParamsPhreeqcLgK)
{
    data["A1"] = obj.A1;
    data["A2"] = obj.A2;
    data["A3"] = obj.A3;
    data["A4"] = obj.A4;
    data["A5"] = obj.A5;
    data["A6"] = obj.A6;
    data["Pr"] = obj.Pr;
}

REAKTORO_DATA_DECODE_DEFINE(ReactionStandardThermoModelParamsPhreeqcLgK)
{
    data.required("A1").to(obj.A1);
    data.required("A2").to(obj.A2);
    data.required("A3").to(obj.A3);
    data.required("A4").to(obj.A4);
    data.required("A5").to(obj.A5);
    data.required("A6").to(obj.A6);
    data.optional("Pr").to(obj.Pr);
}

//----------------------------------------------------------------------

REAKTORO_DATA_ENCODE_DEFINE(ReactionStandardThermoModelParamsVantHoff)
{
    data["lgKr"] = obj.lgKr;
    data["dHr"]  = obj.dHr;
    data["Tr"]   = obj.Tr;
    data["Pr"]   = obj.Pr;
}

REAKTORO_DATA_DECODE_DEFINE(ReactionStandardThermoModelParamsVantHoff)
{
    data.required("lgKr").to(obj.lgKr);
    data.required("dHr").to(obj.dHr);
    data.optional("Tr").to(obj.Tr);
    data.optional("Pr").to(obj.Pr);
}

//======================================================================
// StandardVolumeModelParams Types
//======================================================================

REAKTORO_DATA_ENCODE_DEFINE(StandardVolumeModelParamsConstant)
{
    data["V0"] = obj.V0;
}

REAKTORO_DATA_DECODE_DEFINE(StandardVolumeModelParamsConstant)
{
    data.required("V0").to(obj.V0);
}

} // namespace Reaktoro
