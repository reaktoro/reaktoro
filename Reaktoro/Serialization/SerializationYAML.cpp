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

// Reaktoro includes (StandardThermoModels)
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
// StandardThermoModel
//=====================================================================================================================

auto operator<<(yaml& node, const StandardThermoModelParamsHKF& obj) -> void
{
    auto child = node["HKF"]; // use shorter name instead of full type name (just the last distinct words)!
    child["Gf"] = obj.Gf;
    child["Hf"] = obj.Hf;
    child["Sr"] = obj.Sr;
    child["a1"] = obj.a1;
    child["a2"] = obj.a2;
    child["a3"] = obj.a3;
    child["a4"] = obj.a4;
    child["c1"] = obj.c1;
    child["c2"] = obj.c2;
    child["wref"] = obj.wref;
    child["charge"] = obj.charge;
    child["formula"] = obj.formula;
}

auto operator>>(const yaml& node, StandardThermoModelParamsHKF& obj) -> void
{
    auto child = node.at("HKF"); // use shorter name instead of full type name (just the last distinct words)!
    child.at("Gf").to(obj.Gf);
    child.at("Hf").to(obj.Hf);
    child.at("Sr").to(obj.Sr);
    child.at("a1").to(obj.a1);
    child.at("a2").to(obj.a2);
    child.at("a3").to(obj.a3);
    child.at("a4").to(obj.a4);
    child.at("c1").to(obj.c1);
    child.at("c2").to(obj.c2);
    child.at("wref").to(obj.wref);
    child.at("charge").to(obj.charge);
    child.at("formula").to(obj.formula);
}

//=====================================================================================================================

auto operator<<(yaml& node, const StandardThermoModelParamsHollandPowell& obj) -> void
{
    auto child = node["HollandPowell"]; // use shorter name instead of full type name (just the last distinct words)!
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
    auto child = node.at("HollandPowell"); // use shorter name instead of full type name (just the last distinct words)!
    child.at("Gf").to(obj.Gf);
    child.at("Hf").to(obj.Hf);
    child.at("Sr").to(obj.Sr);
    child.at("Vr").to(obj.Vr);
    child.at("a").to(obj.a);
    child.at("b").to(obj.b);
    child.at("c").to(obj.c);
    child.at("d").to(obj.d);
    child.at("alpha0").to(obj.alpha0);
    child.at("kappa0").to(obj.kappa0);
    child.at("kappa0p").to(obj.kappa0p);
    child.at("kappa0pp").to(obj.kappa0pp);
    child.at("numatoms").to(obj.numatoms);
    child.at("Tcr").to(obj.Tcr);
    child.at("Smax").to(obj.Smax);
    child.at("Vmax").to(obj.Vmax);
    child.at("Tmax").to(obj.Tmax);
}

//=====================================================================================================================

auto operator<<(yaml& node, const StandardThermoModelParamsMaierKelley& obj) -> void
{
    auto child = node["MaierKelley"]; // use shorter name instead of full type name (just the last distinct words)!
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
    auto child = node.at("MaierKelley"); // use shorter name instead of full type name (just the last distinct words)!
    child.at("Gf").to(obj.Gf);
    child.at("Hf").to(obj.Hf);
    child.at("Sr").to(obj.Sr);
    child.at("Vr").to(obj.Vr);
    child.at("a").to(obj.a);
    child.at("b").to(obj.b);
    child.at("c").to(obj.c);
    child.at("Tmax").to(obj.Tmax);
}

//=====================================================================================================================

auto operator<<(yaml& node, const StandardThermoModelParamsMineralHKF& obj) -> void
{
    auto child = node["MineralHKF"]; // use shorter name instead of full type name (just the last distinct words)!
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
    auto child = node.at("MineralHKF"); // use shorter name instead of full type name (just the last distinct words)!
    child.at("Gf").to(obj.Gf);
    child.at("Hf").to(obj.Hf);
    child.at("Sr").to(obj.Sr);
    child.at("Vr").to(obj.Vr);
    child.at("ntr").to(obj.ntr);
    child.at("a").to(obj.a);
    child.at("b").to(obj.b);
    child.at("c").to(obj.c);
    child.at("Ttr").to(obj.Ttr);
    child.at("Htr").to(obj.Htr);
    child.at("Vtr").to(obj.Vtr);
    child.at("dPdTtr").to(obj.dPdTtr);
    child.at("Tmax").to(obj.Tmax);
}

//=====================================================================================================================

auto operator<<(yaml& node, const StandardThermoModelParamsWaterHKF& obj) -> void
{
    auto child = node["WaterHKF"]; // use shorter name instead of full type name (just the last distinct words)!
    node["Ttr"] = obj.Ttr;
    node["Str"] = obj.Str;
    node["Gtr"] = obj.Gtr;
    node["Htr"] = obj.Htr;
}

auto operator>>(const yaml& node, StandardThermoModelParamsWaterHKF& obj) -> void
{
    auto child = node.at("WaterHKF"); // use shorter name instead of full type name (just the last distinct words)!
    child.at("Ttr").to(obj.Ttr);
    child.at("Str").to(obj.Str);
    child.at("Gtr").to(obj.Gtr);
    child.at("Htr").to(obj.Htr);
}

} // namespace YAML
