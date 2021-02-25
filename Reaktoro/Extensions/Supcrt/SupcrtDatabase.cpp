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

#include "SupcrtDatabase.hpp"

// CMakeRC includes
#include <cmrc/cmrc.hpp>

CMRC_DECLARE(ReaktoroDatabases);

// C++ includes
#include <clocale>
#include <sstream>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/AggregateState.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Extensions/Supcrt/SupcrtParams.hpp>
#include <Reaktoro/Extensions/Supcrt/SupcrtEngine.hpp>

// pugixml includes
#include <pugixml.hpp>
using namespace pugi;

namespace Reaktoro {
namespace {

/// The calculator of standard thermodynamic properties based on SUPCRT
SupcrtEngine engine; // TODO: This should be a data member in SupcrtDatabase.

/// Return the standard thermodynamic property function of solvent water using the HKF model.
auto createStandardThermoModel(const SupcrtParamsAqueousSolventHKF& params) -> StandardThermoModel
{
    return [=](real T, real P) { return engine.props(T, P, params); };
}

/// Return the standard thermodynamic property function of an aqueous solute using the HKF model.
auto createStandardThermoModel(const SupcrtParamsAqueousSoluteHKF& params) -> StandardThermoModel
{
    return [=](real T, real P) { return engine.props(T, P, params); };
}

/// Return the standard thermodynamic property function of a fluid species using the Maier-Kelly model.
auto createStandardThermoModel(const SupcrtParamsMaierKelly& params) -> StandardThermoModel
{
    return [=](real T, real P) { return engine.props(T, P, params); };
}

/// Return the standard thermodynamic property function of a mineral species using the Maier-Kelly-HKF model.
auto createStandardThermoModel(const SupcrtParamsMaierKellyHKF& params) -> StandardThermoModel
{
    return [=](real T, real P) { return engine.props(T, P, params); };
}

/// Return the interger value stored in a xml node if not empty.
auto as_int(const xml_node& node, const char* childname, int if_empty=999999) -> int
{
    if(node.child(childname).text().empty())
        return if_empty;
    return node.child(childname).text().as_int();
}

/// Return the double value stored in a xml node if not empty.
auto as_double(const xml_node& node, const char* childname, double if_empty=999999) -> double
{
    if(node.child(childname).text().empty())
        return if_empty;
    return node.child(childname).text().as_double();
}

/// Return the text value stored in a xml node if not empty.
auto as_text(const xml_node& node, const char* childname, String if_empty="999999") -> String
{
    if(node.child(childname).text().empty())
        return if_empty;
    return node.child(childname).text().get();
}

/// Return the text value stored in a xml node if not empty (in lower case).
auto as_lowercase_text(const xml_node& node, const char* childname, String if_empty="999999") -> String
{
    return lowercase(as_text(node, childname, if_empty));
}

/// Return the HKF parameters for the aqueous solvent species stored in an xml node.
auto parseParamsAqueousSolventHKF(const xml_node& node) -> SupcrtParamsAqueousSolventHKF
{
    return {};
}

/// Return the HKF parameters for an aqueous solute stored in an xml node.
auto parseParamsAqueousSoluteHKF(const xml_node& node) -> SupcrtParamsAqueousSoluteHKF
{
    SupcrtParamsAqueousSoluteHKF data;

    data.name = as_text(node, "Name");

    data.charge = as_double(node, "Charge");

    const auto hkf = node.child("Thermo").child("HKF");

    data.Gf   = as_double(hkf, "Gf");
    data.Hf   = as_double(hkf, "Hf");
    data.Sr   = as_double(hkf, "Sr");
    data.a1   = as_double(hkf, "a1");
    data.a2   = as_double(hkf, "a2");
    data.a3   = as_double(hkf, "a3");
    data.a4   = as_double(hkf, "a4");
    data.c1   = as_double(hkf, "c1");
    data.c2   = as_double(hkf, "c2");
    data.wref = as_double(hkf, "wref");

    return data;
}

/// Return the Maier-Kelly parameters for a fluid species stored in an xml node.
auto parseParamsMaierKelly(const xml_node& node) -> SupcrtParamsMaierKelly
{
    SupcrtParamsMaierKelly data;

    data.name = as_text(node, "Name");

    const auto hkf = node.child("Thermo").child("HKF");

    data.Gf   = as_double(hkf, "Gf");
    data.Hf   = as_double(hkf, "Hf");
    data.Sr   = as_double(hkf, "Sr");
    data.a    = as_double(hkf, "a");
    data.b    = as_double(hkf, "b");
    data.c    = as_double(hkf, "c");
    data.Tmax = as_double(hkf, "Tmax");

    return data;
}

/// Return the Maier-Kelly-HKF parameters for a mineral species stored in an xml node.
auto parseParamsMaierKellyHKF(const xml_node& node) -> SupcrtParamsMaierKellyHKF
{
    SupcrtParamsMaierKellyHKF data;

    data.name = as_text(node, "Name");

    const auto hkf = node.child("Thermo").child("HKF");

    data.Gf      = as_double(hkf, "Gf");
    data.Hf      = as_double(hkf, "Hf");
    data.Sr      = as_double(hkf, "Sr");
    data.Vr      = as_double(hkf, "Vr");
    data.nptrans = as_int(hkf, "NumPhaseTrans");
    data.Tmax    = as_double(hkf, "Tmax");

    if(data.nptrans == 0)
    {
        data.a.push_back(as_double(hkf, "a"));
        data.b.push_back(as_double(hkf, "b"));
        data.c.push_back(as_double(hkf, "c"));
    }
    else
    {
        for(int i = 0; i <= data.nptrans; ++i)
        {
            std::stringstream str;
            str << "TemperatureRange" << i;

            auto temperature_range = hkf.child(str.str().c_str());

            data.a.push_back(as_double(temperature_range, "a"));
            data.b.push_back(as_double(temperature_range, "b"));
            data.c.push_back(as_double(temperature_range, "c"));

            if(i < data.nptrans)
            {
                data.Ttr.push_back(as_double(temperature_range, "Ttr"));

                // Set zero the non-available transition values
                data.Htr.push_back(as_double(temperature_range, "Htr", 0.0));
                data.Vtr.push_back(as_double(temperature_range, "Vtr", 0.0));
                data.dPdTtr.push_back(as_double(temperature_range, "dPdTtr", 0.0));
            }
        }
    }

    return data;
}

/// Return the SUPCRT parameters of the species in the xml node as a std::any object.
auto parseParams(const xml_node& node) -> std::any
{
    const auto formula = as_text(node, "Formula");
    const auto type = as_lowercase_text(node, "Type");
    if(type == "aqueous" && formula == "H2O") return parseParamsAqueousSolventHKF(node);
    if(type == "aqueous") return parseParamsAqueousSoluteHKF(node);
    if(type == "gaseous") return parseParamsMaierKelly(node);
    if(type == "mineral") return parseParamsMaierKellyHKF(node);
    return parseParamsMaierKelly(node); // assume Maier-Kelly data by default
}

/// Return the pairs of element symbols and their coefficients in an xml node.
auto parseElements(const xml_node& node) -> Map<String, double>
{
    const auto nodetext = as_text(node, "Elements");
    auto words = split(nodetext, "()");
    Map<String, double> symbols;
    for(unsigned i = 0; i < words.size(); i += 2)
        symbols.emplace(words[i], tofloat(words[i + 1]));
    return symbols;
}

/// Return the aggregate state of the species in the xml node.
auto parseAggregateState(const xml_node& node)
{
    const auto type = as_lowercase_text(node, "Type");
    if(type == "aqueous") return AggregateState::Aqueous;
    if(type == "gaseous") return AggregateState::Gas;
    if(type == "liquid") return AggregateState::Liquid;
    if(type == "mineral") return AggregateState::Solid;
    return AggregateState::Undefined;
}

/// Return the standard thermo props function of the species node based on SUPCRT.
auto parseStandardThermoModel(const xml_node& node) -> StandardThermoModel
{
    const auto formula = as_text(node, "Formula");
    const auto type = as_lowercase_text(node, "Type");
    if(type == "aqueous" && formula == "H2O") return createStandardThermoModel(parseParamsAqueousSolventHKF(node));
    if(type == "aqueous") return createStandardThermoModel(parseParamsAqueousSoluteHKF(node));
    if(type == "gaseous") return createStandardThermoModel(parseParamsMaierKelly(node));
    if(type == "mineral") return createStandardThermoModel(parseParamsMaierKellyHKF(node));
    return createStandardThermoModel(parseParamsMaierKelly(node)); // assume Maier-Kelly data by default
}

auto parseSpecies(const xml_node& node) -> Species
{
    Species species;
    species = species.withName(as_text(node, "Name"));
    species = species.withFormula(as_text(node, "Formula"));
    species = species.withElements(parseElements(node));
    species = species.withCharge(as_double(node, "Charge", 0.0));
    species = species.withAggregateState(parseAggregateState(node));
    species = species.withStandardThermoModel(parseStandardThermoModel(node));
    species = species.withAttachedData(parseParams(node));
    return species;
}

auto parseDatabase(const xml_document& doc) -> SupcrtDatabase
{
    xml_node database = doc.child("Database");

    SupcrtDatabase db;
    for(xml_node node : doc.child("Database").children("Species"))
        db.addSpecies(parseSpecies(node));

    return db;
}

/// Return the contents of the embedded SUPCRT database with given name (or empty)
auto getSupcrtDatabaseContent(String name) -> String
{
    error(!oneof(name,
        "supcrt98.xml",
        "supcrt07.xml",
        "supcrt98-organics.xml",
        "supcrt07-organics.xml"),
        "Could not load embedded Supcrt database file with name `", name, "`. ",
        "The currently supported names are: \n"
        "    - supcrt98.xml          \n",
        "    - supcrt07.xml          \n",
        "    - supcrt98-organics.xml \n",
        "    - supcrt07-organics.xml \n",
        "");
    auto fs = cmrc::ReaktoroDatabases::get_filesystem();
    auto contents = fs.open("databases/supcrt/" + name);
    return String(contents.begin(), contents.end());
}

} // namespace

/// An auxiliary type to change locale and ensure its return to original.
/// This is needed to avoid certain issues with pugixml related to how decimal numbers are represented in different languages.
struct ChangeLocale
{
    const String old_locale;

    explicit ChangeLocale(const char* new_locale) : old_locale(std::setlocale(LC_NUMERIC, nullptr))
    {
        std::setlocale(LC_NUMERIC, new_locale);
    }

    ~ChangeLocale()
    {
        std::setlocale(LC_NUMERIC, old_locale.c_str());
    }
};

SupcrtDatabase::SupcrtDatabase()
{}

SupcrtDatabase::SupcrtDatabase(String name)
: SupcrtDatabase(SupcrtDatabase::withName(name))
{}

auto SupcrtDatabase::withName(String name) ->  SupcrtDatabase
{
    // Change locale to C at construction and reset at destruction.
    const auto guard = ChangeLocale("C");

    // Get the text content of an embedded SUPCRT database
    String text = getSupcrtDatabaseContent(name);

    // Assert the given database name is valid
    Assert(text.size(), "Could not create a SupcrtDatabase object.",
        "There is no embedded SUPCRT database with name " + name);

    // Parse the embedded database text content
    xml_document doc;
    auto result = doc.load_string(text.c_str());

    // Assert the xml parsing is successfull
    error(!result, "Could not create a SupcrtDatabase object. ",
        "There was an error parsing the embedded SUPCRT database with name " + name + ". "
        "The parsing error was:\n" + String(result.description()));

    // Parse the xml document into a SupcrtDatabase object
    return parseDatabase(doc);
}

auto SupcrtDatabase::fromFile(String path) ->  SupcrtDatabase
{
    // Change locale to C at construction and reset at destruction.
    const auto guard = ChangeLocale("C");

    // Parse the xml database file
    xml_document doc;
    auto result = doc.load_file(path.c_str());

    // Assert the xml parsing is successfull
    Assert(result, "Could not create a SupcrtDatabase object.",
        "There was an error parsing the SUPCRT database at given path: " + path + ". The parsing error was:\n" + String(result.description()));

    // Parse the xml document into a SupcrtDatabase object
    return parseDatabase(doc);
}

} // namespace Reaktoro

