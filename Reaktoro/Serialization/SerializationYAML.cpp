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

// Reaktoro includes (StandardThermoModel)
#include <Reaktoro/Thermodynamics/Standard/StandardThermoModelHKF.hpp>
#include <Reaktoro/Thermodynamics/Standard/StandardThermoModelHollandPowell.hpp>
#include <Reaktoro/Thermodynamics/Standard/StandardThermoModelMaierKelley.hpp>
#include <Reaktoro/Thermodynamics/Standard/StandardThermoModelMineralHKF.hpp>
#include <Reaktoro/Thermodynamics/Standard/StandardThermoModelWaterHKF.hpp>

namespace YAML {

//======================================================================
// Common
//======================================================================

// Needed to avoid ambiguous overloaded operator '<<' error.
auto operator<<(Node& node, double obj) -> void
{
    node = obj;
}

// Needed to avoid ambiguous overloaded operator '>>' error.
auto operator>>(const Node& node, double& obj) -> void
{
    obj = node.as<double>();
}

auto operator<<(Node& node, const real& obj) -> void
{
    node = obj.val();
}

auto operator>>(const Node& node, real& obj) -> void
{
    obj = node.as<double>();
}

//=====================================================================================================================
// Core
//=====================================================================================================================

auto operator<<(Node& node, const ChemicalFormula& obj) -> void
{
    node = obj.str();
}

auto operator>>(const Node& node, ChemicalFormula& obj) -> void
{
    obj = ChemicalFormula(node.as<std::string>());
}

//=====================================================================================================================

auto operator<<(Node& node, const ChemicalSystem& obj) -> void
{
    errorif(true, "`auto operator<<(Node& node, const ChemicalSystem& obj) -> void` not implemented!");
}

auto operator>>(const Node& node, ChemicalSystem& obj) -> void
{
    errorif(true,  "`auto operator<<=(ChemicalSystem& obj) -> void` not implemented!");
}

//=====================================================================================================================

auto operator<<(Node& node, const Element& obj) -> void
{
    appendIfNotDefault(node, "Symbol",            obj.symbol());
    appendIfNotDefault(node, "Name",              obj.name());
    appendIfNotDefault(node, "AtomicNumber",      obj.atomicNumber());
    appendIfNotDefault(node, "AtomicWeight",      obj.atomicWeight());
    appendIfNotDefault(node, "MolarMass",         obj.molarMass());
    appendIfNotDefault(node, "Electronegativity", obj.electronegativity());
    appendIfNotDefault(node, "Tags",              obj.tags());
}

auto operator>>(const Node& node, Element& obj) -> void
{
    Element::Args args;
    args.symbol            = get(node, "Symbol", String{});
    args.name              = get(node, "Name", String{});
    args.atomic_number     = get(node, "AtomicNumber", 0);
    args.atomic_weight     = get(node, "AtomicWeight", 0.0);
    args.electronegativity = get(node, "Electronegativity", 0.0);
    args.tags              = get(node, "Tags", Strings{});

    obj = Element(args);
}

//=====================================================================================================================

auto operator<<(Node& node, const Param& obj) -> void
{
    node = obj.value();
}

auto operator>>(const Node& node, Param& obj) -> void
{
    obj = node.as<double>();
}

//=====================================================================================================================

auto operator<<(Node& node, const Params& obj) -> void
{
    node = obj.data();
}

auto operator>>(const Node& node, Params& obj) -> void
{
    auto values = node.as<Vec<double>>();
    obj = Params(values.begin(), values.end());
}

//=====================================================================================================================

auto operator<<(Node& node, const Phase& obj) -> void
{
    errorif(true, "`auto operator<<(Node& node, const Phase& obj) -> void` not implemented!");
}

auto operator>>(const Node& node, Phase& obj) -> void
{
    errorif(true,  "`auto operator<<=(Phase& obj) -> void` not implemented!");
}

//=====================================================================================================================

auto operator<<(Node& node, const Species& obj) -> void
{

}

auto operator>>(const Node& node, Species& obj) -> void
{

}

//=====================================================================================================================
// StandardThermoModel
//=====================================================================================================================

auto operator<<(Node& node, const StandardThermoModelParamsHKF& obj) -> void
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

auto operator>>(const Node& node, StandardThermoModelParamsHKF& obj) -> void
{
    YAML::set(node, "Gf", obj.Gf);
    YAML::set(node, "Hf", obj.Hf);
    YAML::set(node, "Sr", obj.Sr);
    YAML::set(node, "a1", obj.a1);
    YAML::set(node, "a2", obj.a2);
    YAML::set(node, "a3", obj.a3);
    YAML::set(node, "a4", obj.a4);
    YAML::set(node, "c1", obj.c1);
    YAML::set(node, "c2", obj.c2);
    YAML::set(node, "wref", obj.wref);
    YAML::set(node, "charge", obj.charge);
    YAML::set(node, "formula", obj.formula);
}

//=====================================================================================================================

auto operator<<(Node& node, const StandardThermoModelParamsHollandPowell& obj) -> void
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

auto operator>>(const Node& node, StandardThermoModelParamsHollandPowell& obj) -> void
{
    YAML::set(node, "Gf",       obj.Gf);
    YAML::set(node, "Hf",       obj.Hf);
    YAML::set(node, "Sr",       obj.Sr);
    YAML::set(node, "Vr",       obj.Vr);
    YAML::set(node, "a",        obj.a);
    YAML::set(node, "b",        obj.b);
    YAML::set(node, "c",        obj.c);
    YAML::set(node, "d",        obj.d);
    YAML::set(node, "alpha0",   obj.alpha0);
    YAML::set(node, "kappa0",   obj.kappa0);
    YAML::set(node, "kappa0p",  obj.kappa0p);
    YAML::set(node, "kappa0pp", obj.kappa0pp);
    YAML::set(node, "numatoms", obj.numatoms);
    YAML::set(node, "Tcr",      obj.Tcr);
    YAML::set(node, "Smax",     obj.Smax);
    YAML::set(node, "Vmax",     obj.Vmax);
    YAML::set(node, "Tmax",     obj.Tmax);
}

//=====================================================================================================================

auto operator<<(Node& node, const StandardThermoModelParamsMaierKelley& obj) -> void
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

auto operator>>(const Node& node, StandardThermoModelParamsMaierKelley& obj) -> void
{
    YAML::set(node, "Gf",   obj.Gf);
    YAML::set(node, "Hf",   obj.Hf);
    YAML::set(node, "Sr",   obj.Sr);
    YAML::set(node, "Vr",   obj.Vr);
    YAML::set(node, "a",    obj.a);
    YAML::set(node, "b",    obj.b);
    YAML::set(node, "c",    obj.c);
    YAML::set(node, "Tmax", obj.Tmax);
}

//=====================================================================================================================

auto operator<<(Node& node, const StandardThermoModelParamsMineralHKF& obj) -> void
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

auto operator>>(const Node& node, StandardThermoModelParamsMineralHKF& obj) -> void
{
    YAML::set(node, "Gf",     obj.Gf);
    YAML::set(node, "Hf",     obj.Hf);
    YAML::set(node, "Sr",     obj.Sr);
    YAML::set(node, "Vr",     obj.Vr);
    YAML::set(node, "ntr",    obj.ntr);
    YAML::set(node, "a",      obj.a);
    YAML::set(node, "b",      obj.b);
    YAML::set(node, "c",      obj.c);
    YAML::set(node, "Ttr",    obj.Ttr);
    YAML::set(node, "Htr",    obj.Htr);
    YAML::set(node, "Vtr",    obj.Vtr);
    YAML::set(node, "dPdTtr", obj.dPdTtr);
    YAML::set(node, "Tmax",   obj.Tmax);
}

//=====================================================================================================================

auto operator<<(Node& node, const StandardThermoModelParamsWaterHKF& obj) -> void
{
    node["Ttr"] = obj.Ttr;
    node["Str"] = obj.Str;
    node["Gtr"] = obj.Gtr;
    node["Htr"] = obj.Htr;
}

auto operator>>(const Node& node, StandardThermoModelParamsWaterHKF& obj) -> void
{
    YAML::set(node, "Ttr", obj.Ttr);
    YAML::set(node, "Str", obj.Str);
    YAML::set(node, "Gtr", obj.Gtr);
    YAML::set(node, "Htr", obj.Htr);
}

} // namespace YAML
