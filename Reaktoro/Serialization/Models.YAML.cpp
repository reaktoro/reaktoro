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

#include "Models.YAML.hpp"

// Reaktoro includes
#include <Reaktoro/Models/ReactionThermoModelConstLgK.hpp>
#include <Reaktoro/Models/ReactionThermoModelGemsLgK.hpp>
#include <Reaktoro/Models/ReactionThermoModelPhreeqcLgK.hpp>
#include <Reaktoro/Models/ReactionThermoModelVantHoff.hpp>
#include <Reaktoro/Models/StandardThermoModelHKF.hpp>
#include <Reaktoro/Models/StandardThermoModelHollandPowell.hpp>
#include <Reaktoro/Models/StandardThermoModelMaierKelley.hpp>
#include <Reaktoro/Models/StandardThermoModelMineralHKF.hpp>
#include <Reaktoro/Models/StandardThermoModelWaterHKF.hpp>
#include <Reaktoro/Serialization/Common.YAML.hpp>
#include <Reaktoro/Serialization/Core.YAML.hpp>

namespace Reaktoro {

REAKTORO_YAML_ENCODE_DEFINE(ReactionThermoModelParamsConstLgK)
{
    node["lgKr"] = obj.lgKr;
    node["Pr"]   = obj.Pr;
}

REAKTORO_YAML_DECODE_DEFINE(ReactionThermoModelParamsConstLgK)
{
    node.at("lgKr").to(obj.lgKr);
    node.at("Pr").to(obj.Pr);
}

//=====================================================================================================================

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
    node.at("A0").to(obj.A0);
    node.at("A1").to(obj.A1);
    node.at("A2").to(obj.A2);
    node.at("A3").to(obj.A3);
    node.at("A4").to(obj.A4);
    node.at("A5").to(obj.A5);
    node.at("A6").to(obj.A6);
    node.at("Pr").to(obj.Pr);
}

//=====================================================================================================================

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
    node.at("A1").to(obj.A1);
    node.at("A2").to(obj.A2);
    node.at("A3").to(obj.A3);
    node.at("A4").to(obj.A4);
    node.at("A5").to(obj.A5);
    node.at("A6").to(obj.A6);
    node.at("Pr").to(obj.Pr);
}

//=====================================================================================================================

REAKTORO_YAML_ENCODE_DEFINE(ReactionThermoModelParamsVantHoff)
{
    node["lgKr"] = obj.lgKr;
    node["dHr"]  = obj.dHr;
    node["Tr"]   = obj.Tr;
    node["Pr"]   = obj.Pr;
}

REAKTORO_YAML_DECODE_DEFINE(ReactionThermoModelParamsVantHoff)
{
    node.at("lgKr").to(obj.lgKr);
    node.at("dHr").to(obj.dHr);
    node.at("Tr").to(obj.Tr);
    node.at("Pr").to(obj.Pr);
}

//=====================================================================================================================

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

//=====================================================================================================================

REAKTORO_YAML_ENCODE_DEFINE(StandardThermoModelParamsHollandPowell)
{
    node["Gf"]       = obj.Gf;
    node["Hf"]       = obj.Hf;
    node["Sr"]       = obj.Sr;
    node["Vr"]       = obj.Vr;
    node["a"]        = obj.a;
    node["b"]        = obj.b;
    node["c"]        = obj.c;
    node["d"]        = obj.d;
    node["alpha0"]   = obj.alpha0;
    node["kappa0"]   = obj.kappa0;
    node["kappa0p"]  = obj.kappa0p;
    node["kappa0pp"] = obj.kappa0pp;
    node["numatoms"] = obj.numatoms;
    node["Tmax"]     = obj.Tmax;
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
    node.at("alpha0").to(obj.alpha0);
    node.at("kappa0").to(obj.kappa0);
    node.at("kappa0p").to(obj.kappa0p);
    node.at("kappa0pp").to(obj.kappa0pp);
    node.at("numatoms").to(obj.numatoms);
    node.at("Tmax").to(obj.Tmax);
}

//=====================================================================================================================

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

//=====================================================================================================================

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

//=====================================================================================================================

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

} // namespace YAML
