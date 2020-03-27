// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include "DatabaseSupcrt.hpp"

// C++ includes
#include <clocale>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Databases/EmbeddedDatabases.hpp>
#include <Reaktoro/Thermodynamics/Models/SpeciesThermoStateHKF.hpp>

// miniz includes
#include <miniz/zip_file.hpp>

// pugixml includes
#include <pugixml.hpp>
using namespace pugi;

namespace Reaktoro {
namespace {

auto errorNonExistentSpecies(std::string type, std::string name) -> void
{
    Exception exception;
    exception.error << "Cannot get an instance of the " << type << " species `" << name << "` in the database.";
    exception.reason << "There is no such species in the database.";
    RaiseError(exception);
}

auto parseDissociation(std::string dissociation) -> std::map<std::string, double>
{
    std::map<std::string, double> equation;
    auto words = split(dissociation, " ");
    for(const auto& word : words)
    {
        auto pair = split(word, ":");
        equation.emplace(pair[1], tofloat(pair[0]));
    }
    return equation;
}

auto as_int(const xml_node& node, const char* childname, int if_empty=999999) -> int
{
    if(node.child(childname).text().empty())
        return if_empty;
    return node.child(childname).text().as_int();
}

auto as_double(const xml_node& node, const char* childname, double if_empty=999999) -> double
{
    if(node.child(childname).text().empty())
        return if_empty;
    return node.child(childname).text().as_double();
}

auto as_text(const xml_node& node, const char* childname, std::string if_empty="999999") -> std::string
{
    if(node.child(childname).text().empty())
        return if_empty;
    return node.child(childname).text().get();
}

auto as_lowercase_text(const xml_node& node, const char* childname, std::string if_empty="999999") -> std::string
{
    return lowercase(as_text(node, childname, if_empty));
}

auto parseElement(const xml_node& node) -> Element
{
    Element element;
    element = element.withName(as_text(node, "Name"));
    element = element.withMolarMass(as_double(node, "MolarMass") * 1.0e-3); // convert from g/mol to kg/mol
    return element;
}

auto parseParamsAqueousSoluteHKF(const xml_node& node) -> ParamsAqueousSoluteHKF
{
    ParamsAqueousSoluteHKF data;

    data.name = as_text(node, "Name");

    data.charge = as_double(node, "Charge");
    data.dissociation = parseDissociation(as_text(node, "Dissociation"));

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

auto parseParamsMaierKelly(const xml_node& node) -> ParamsMaierKelly
{
    ParamsMaierKelly data;

    data.name = as_text(node, "Name");

    data.Tcr = as_double(node, "CriticalTemperature");
    data.Pcr = as_double(node, "CriticalPressure") * barToPascal; // convert from bar to Pa
    data.acentric_factor = as_double(node, "AcentricFactor");

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

auto parseParamsMaierKellyHKF(const xml_node& node) -> ParamsMaierKellyHKF
{
    ParamsMaierKellyHKF data;

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

auto parseElementalFormula(const xml_node& node, const ElementMap& element_map) -> std::map<Element, double>
{
    std::string formula = as_text(node, "Elements");
    std::map<Element, double> elements;
    auto words = split(formula, "()");
    for(unsigned i = 0; i < words.size(); i += 2)
    {
        Assert(element_map.count(words[i]),
            "Cannot parse the elemental formula `" + formula + "`.",
            "The element `" + words[i] + "` is not in the database.");
        elements.emplace(element_map.at(words[i]), tofloat(words[i + 1]));
    }
    if(!node.child("Charge").empty())
    {
        double charge = node.child("Charge").text().as_double();
        elements.emplace(element_map.at("Z"), charge);
    }
    return elements;
}

auto parseSpeciesData(const xml_node& node) -> std::any
{
    const auto type = as_lowercase_text(node, "Type");
    if(type == "aqueous") return parseParamsAqueousSoluteHKF(node);
    if(type == "gaseous") return parseParamsMaierKelly(node);
    if(type == "mineral") return parseParamsMaierKellyHKF(node);
    return parseParamsMaierKelly(node); // assume Maier-Kelly data by default
}

auto parseSpecies(const xml_node& node, const ElementMap& element_map) -> Species
{
    Species species;
    species = species.withName(as_text(node, "Name"));
    species = species.withFormula(as_text(node, "Formula"));
    species = species.withElements(parseElementalFormula(node, element_map));
    species = species.withType(as_lowercase_text(node, "Type"));
    species = species.withData(parseSpeciesData(node));
    return species;
}

auto parseElementMap(const xml_document& doc) -> ElementMap
{
    ElementMap element_map;
    for(xml_node node : doc.child("Database").children("Element"))
    {
        Element element = parseElement(node);
        element_map[element.name()] = element;
    }
    element_map["Z"] = Element().withName("Z");
    return element_map;
}

auto parseSpeciesMap(const xml_document& doc, const ElementMap& element_map) -> SpeciesMap
{
    SpeciesMap species_map;
    for(xml_node node : doc.child("Database").children("Species"))
    {
        Species species = parseSpecies(node, element_map);
        species_map[species.name()] = species;
    }
    return species_map;
}

auto parseDatabase(const xml_document& doc) -> DatabaseSupcrt
{
    xml_node database = doc.child("Database");

    ElementMap element_map = parseElementMap(doc);
    SpeciesMap species_map = parseSpeciesMap(doc, element_map);

    DatabaseSupcrt db;
    db.setElements(element_map);
    db.setSpecies(species_map);
    return db;
}

} // namespace

/// An auxiliary type to change locale and ensure its return to original.
/// This is needed to avoid certain issues with pugixml related to how decimal numbers are represented in different languages.
struct ChangeLocale
{
    const std::string old_locale;

    explicit ChangeLocale(const char* new_locale) : old_locale(std::setlocale(LC_NUMERIC, nullptr))
    {
        std::setlocale(LC_NUMERIC, new_locale);
    }

    ~ChangeLocale()
    {
        std::setlocale(LC_NUMERIC, old_locale.c_str());
    }
};

DatabaseSupcrt::DatabaseSupcrt()
{}

auto DatabaseSupcrt::withName(std::string name) ->  DatabaseSupcrt
{
    // Change locale to C at construction and reset at destruction.
    const auto guard = ChangeLocale("C");

    // Get the text content of a SUPCRT database
    std::string text = supcrtEmbeddedDatabaseTextContent(name);

    // Assert the given database name is valid
    Assert(text.size(), "Could not create a DatabaseSupcrt object.",
        "There is no embedded SUPCRT database with name " + name);

    // Parse the embedded database text content
    xml_document doc;
    auto result = doc.load_string(text.c_str());

    // Assert the xml parsing is successfull
    Assert(result, "Could not create a DatabaseSupcrt object.",
        "There was an error parsing the embedded SUPCRT database with name " + name + ". The parsing error was:\n" + std::string(result.description()));

    // Parse the xml document into a DatabaseSupcrt object
    return parseDatabase(doc);
}

auto DatabaseSupcrt::fromFile(std::string path) ->  DatabaseSupcrt
{
    // Change locale to C at construction and reset at destruction.
    const auto guard = ChangeLocale("C");

    // Parse the xml database file
    xml_document doc;
    auto result = doc.load_file(path.c_str());

    // Assert the xml parsing is successfull
    Assert(result, "Could not create a DatabaseSupcrt object.",
        "There was an error parsing the SUPCRT database at given path: " + path + ". The parsing error was:\n" + std::string(result.description()));

    // Parse the xml document into a DatabaseSupcrt object
    return parseDatabase(doc);
}

} // namespace Reaktoro

