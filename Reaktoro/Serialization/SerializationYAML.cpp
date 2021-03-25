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

#include "SerializationYAML.hpp"

// Reaktoro includes (Common)
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/YAML.hpp>

// Reaktoro includes (Core)
#include <Reaktoro/Core/ChemicalFormula.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Param.hpp>
#include <Reaktoro/Core/Params.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/Species.hpp>

// Reaktoro includes (Models)
#include <Reaktoro/Models/ReactionThermoModelConstLgK.hpp>
#include <Reaktoro/Models/ReactionThermoModelGemsLgK.hpp>
#include <Reaktoro/Models/ReactionThermoModelPhreeqcLgK.hpp>
#include <Reaktoro/Models/ReactionThermoModelVantHoff.hpp>
#include <Reaktoro/Models/StandardThermoModelHKF.hpp>
#include <Reaktoro/Models/StandardThermoModelHollandPowell.hpp>
#include <Reaktoro/Models/StandardThermoModelMaierKelley.hpp>
#include <Reaktoro/Models/StandardThermoModelMineralHKF.hpp>
#include <Reaktoro/Models/StandardThermoModelWaterHKF.hpp>

namespace Reaktoro {

//======================================================================
// Common
//======================================================================

// Needed to avoid ambiguous overloaded operator '<<' error.
auto operator<<(yaml& node, double obj) -> void
{
    node = obj;
}

// Needed to avoid ambiguous overloaded operator '>>' error.
auto operator>>(const yaml& node, double& obj) -> void
{
    obj = node.as<double>();
}

auto operator<<(yaml& node, const real& obj) -> void
{
    node = obj.val();
}

auto operator>>(const yaml& node, real& obj) -> void
{
    obj = node.as<double>();
}

//=====================================================================================================================
// Core
//=====================================================================================================================

auto operator<<(yaml& node, const ChemicalFormula& obj) -> void
{
    node = obj.str();
}

auto operator>>(const yaml& node, ChemicalFormula& obj) -> void
{
    obj = ChemicalFormula(node.as<std::string>());
}

//=====================================================================================================================

auto operator<<(yaml& node, const ChemicalSystem& obj) -> void
{
    errorif(true, "`auto operator<<(yaml& node, const ChemicalSystem& obj) -> void` not implemented!");
}

auto operator>>(const yaml& node, ChemicalSystem& obj) -> void
{
    errorif(true,  "`auto operator<<=(ChemicalSystem& obj) -> void` not implemented!");
}

//=====================================================================================================================

auto operator<<(yaml& node, const Element& obj) -> void
{
    node.appendIfNotDefault("Symbol",            obj.symbol());
    node.appendIfNotDefault("Name",              obj.name());
    node.appendIfNotDefault("AtomicNumber",      obj.atomicNumber());
    node.appendIfNotDefault("AtomicWeight",      obj.atomicWeight());
    node.appendIfNotDefault("MolarMass",         obj.molarMass());
    node.appendIfNotDefault("Electronegativity", obj.electronegativity());
    node.appendIfNotDefault("Tags",              obj.tags());
}

auto operator>>(const yaml& node, Element& obj) -> void
{
    Element::Args args;
    args.symbol            = node["Symbol"].value(String{});
    args.name              = node["Name"].value(String{});
    args.atomic_number     = node["AtomicNumber"].value(0);
    args.atomic_weight     = node["AtomicWeight"].value(0.0);
    args.electronegativity = node["Electronegativity"].value(0.0);
    args.tags              = node["Tags"].value(Strings{});

    obj = Element(args);
}

//=====================================================================================================================

auto operator<<(yaml& node, const Param& obj) -> void
{
    node = obj.value();
}

auto operator>>(const yaml& node, Param& obj) -> void
{
    obj = node.as<double>();
}

//=====================================================================================================================

auto operator<<(yaml& node, const Params& obj) -> void
{
    node = obj.data();
}

auto operator>>(const yaml& node, Params& obj) -> void
{
    auto values = node.as<Vec<double>>();
    obj = Params(values.begin(), values.end());
}

//=====================================================================================================================

auto operator<<(yaml& node, const Phase& obj) -> void
{
    errorif(true, "`auto operator<<(yaml& node, const Phase& obj) -> void` not implemented!");
}

auto operator>>(const yaml& node, Phase& obj) -> void
{
    errorif(true,  "`auto operator<<=(Phase& obj) -> void` not implemented!");
}

//=====================================================================================================================

auto operator<<(yaml& node, const Species& obj) -> void
{

}

auto operator>>(const yaml& node, Species& obj) -> void
{

}

//=====================================================================================================================
// Models
//=====================================================================================================================

auto operator<<(yaml& node, const ReactionThermoModelParamsConstLgK& obj) -> void
{
    node["lgKr"] = obj.lgKr;
    node.appendIfNotDefault("Pr", obj.Pr, 100'000);
}

auto operator>>(const yaml& node, ReactionThermoModelParamsConstLgK& obj) -> void
{
    node.at("lgKr").to(obj.lgKr);
    obj.Pr = node["Pr"].value(100'000);
}

//=====================================================================================================================

auto operator<<(yaml& node, const ReactionThermoModelParamsGemsLgK& obj) -> void
{
    node["A0"] = obj.A0;
    node["A1"] = obj.A1;
    node["A2"] = obj.A2;
    node["A3"] = obj.A3;
    node["A4"] = obj.A4;
    node["A5"] = obj.A5;
    node["A6"] = obj.A6;
    node.appendIfNotDefault("Pr", obj.Pr, 100'000);
}

auto operator>>(const yaml& node, ReactionThermoModelParamsGemsLgK& obj) -> void
{
    node.at("A0").to(obj.A0);
    node.at("A1").to(obj.A1);
    node.at("A2").to(obj.A2);
    node.at("A3").to(obj.A3);
    node.at("A4").to(obj.A4);
    node.at("A5").to(obj.A5);
    node.at("A6").to(obj.A6);
    obj.Pr = node["Pr"].value(100'000);
}

//=====================================================================================================================

auto operator<<(yaml& node, const ReactionThermoModelParamsPhreeqcLgK& obj) -> void
{
    node["A1"] = obj.A1;
    node["A2"] = obj.A2;
    node["A3"] = obj.A3;
    node["A4"] = obj.A4;
    node["A5"] = obj.A5;
    node["A6"] = obj.A6;
    node.appendIfNotDefault("Pr", obj.Pr, 100'000);
}

auto operator>>(const yaml& node, ReactionThermoModelParamsPhreeqcLgK& obj) -> void
{
    node.at("A1").to(obj.A1);
    node.at("A2").to(obj.A2);
    node.at("A3").to(obj.A3);
    node.at("A4").to(obj.A4);
    node.at("A5").to(obj.A5);
    node.at("A6").to(obj.A6);
    obj.Pr = node["Pr"].value(100'000);
}

//=====================================================================================================================

auto operator<<(yaml& node, const ReactionThermoModelParamsVantHoff& obj) -> void
{
    node["lgKr"] = obj.lgKr;
    node["dHr"] = obj.dHr;
    node.appendIfNotDefault("Tr", obj.Tr, 298.15);
    node.appendIfNotDefault("Pr", obj.Pr, 100'000);
}

auto operator>>(const yaml& node, ReactionThermoModelParamsVantHoff& obj) -> void
{
    node.at("lgKr").to(obj.lgKr);
    node.at("dHr").to(obj.dHr);
    obj.Tr = node["Tr"].value(298.15);
    obj.Pr = node["Pr"].value(100'000);
}

//=====================================================================================================================

auto operator<<(yaml& node, const StandardThermoModelParamsHKF& obj) -> void
{
    node["Gf"] = obj.Gf;
    node["Hf"] = obj.Hf;
    node["Sr"] = obj.Sr;
    node["a1"] = obj.a1;
    node["a2"] = obj.a2;
    node["a3"] = obj.a3;
    node["a4"] = obj.a4;
    node["c1"] = obj.c1;
    node["c2"] = obj.c2;
    node["wref"] = obj.wref;
    node["charge"] = obj.charge;
    node["formula"] = obj.formula;
}

auto operator>>(const yaml& node, StandardThermoModelParamsHKF& obj) -> void
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
    node.at("formula").to(obj.formula);
}

//=====================================================================================================================

auto operator<<(yaml& node, const StandardThermoModelParamsHollandPowell& obj) -> void
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
    node["Tcr"]      = obj.Tcr;
    node["Smax"]     = obj.Smax;
    node["Vmax"]     = obj.Vmax;
    node["Tmax"]     = obj.Tmax;
}

auto operator>>(const yaml& node, StandardThermoModelParamsHollandPowell& obj) -> void
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
    node.at("Tcr").to(obj.Tcr);
    node.at("Smax").to(obj.Smax);
    node.at("Vmax").to(obj.Vmax);
    node.at("Tmax").to(obj.Tmax);
}

//=====================================================================================================================

auto operator<<(yaml& node, const StandardThermoModelParamsMaierKelley& obj) -> void
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

auto operator>>(const yaml& node, StandardThermoModelParamsMaierKelley& obj) -> void
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

auto operator<<(yaml& node, const StandardThermoModelParamsMineralHKF& obj) -> void
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

auto operator>>(const yaml& node, StandardThermoModelParamsMineralHKF& obj) -> void
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

auto operator<<(yaml& node, const StandardThermoModelParamsWaterHKF& obj) -> void
{
    node["Ttr"] = obj.Ttr;
    node["Str"] = obj.Str;
    node["Gtr"] = obj.Gtr;
    node["Htr"] = obj.Htr;
}

auto operator>>(const yaml& node, StandardThermoModelParamsWaterHKF& obj) -> void
{
    node.at("Ttr").to(obj.Ttr);
    node.at("Str").to(obj.Str);
    node.at("Gtr").to(obj.Gtr);
    node.at("Htr").to(obj.Htr);
}

} // namespace YAML
