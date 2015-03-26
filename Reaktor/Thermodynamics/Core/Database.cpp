// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "Database.hpp"

// C++ includes
#include <map>
#include <set>
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Common/Constants.hpp>
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Common/Optional.hpp>
#include <Reaktor/Common/SetUtils.hpp>
#include <Reaktor/Common/StringUtils.hpp>
#include <Reaktor/Common/Units.hpp>
#include <Reaktor/Core/Element.hpp>
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktor/Thermodynamics/Species/GaseousSpecies.hpp>
#include <Reaktor/Thermodynamics/Species/MineralSpecies.hpp>

// pugixml includes
#include <Reaktor/pugixml/pugixml.hpp>
using namespace pugi;

namespace Reaktor {
namespace {

/// Auxiliary types for a map of aqueous, gaseous, and mineral species
using ElementMap        = std::map<std::string, Element>;
using AqueousSpeciesMap = std::map<std::string, AqueousSpecies>;
using GaseousSpeciesMap = std::map<std::string, GaseousSpecies>;
using MineralSpeciesMap = std::map<std::string, MineralSpecies>;

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

auto parseReactionThermoProperties(const xml_node& node) -> ReactionThermoProperties
{
    // Get the data values of the children nodes
    std::vector<double> temperatures     = tofloats(node.child("temperatures").text().get());
    std::vector<double> pressures        = tofloats(node.child("pressures").text().get());
    std::vector<double> pk               = tofloats(node.child("pk").text().get());
    std::vector<double> lnk              = tofloats(node.child("lnk").text().get());
    std::vector<double> logk             = tofloats(node.child("logk").text().get());
    std::vector<double> gibbs_energy     = tofloats(node.child("G").text().get());
    std::vector<double> helmholtz_energy = tofloats(node.child("A").text().get());
    std::vector<double> internal_energy  = tofloats(node.child("U").text().get());
    std::vector<double> enthalpy         = tofloats(node.child("H").text().get());
    std::vector<double> entropy          = tofloats(node.child("S").text().get());
    std::vector<double> volume           = tofloats(node.child("V").text().get());
    std::vector<double> heat_capacity    = tofloats(node.child("C").text().get());

    // Convert `pk` to `lnk`, where `pk = -log(k) = -ln(k)/ln(10)`
    const double ln_10 = std::log(10.0);
    if(not pk.empty() and lnk.empty())
        for(unsigned i = 0; i < pk.size(); ++i)
            lnk[i] = -pk[i] * ln_10;

    // Convert `logk` to `lnk`, where `log(k) = ln(k)/ln(10)`
    if(not logk.empty() and lnk.empty())
        for(unsigned i = 0; i < logk.size(); ++i)
            lnk[i] = logk[i] * ln_10;

    // Get the temperature and pressure units
    std::string tunits = node.child("temperatures").attribute("units").as_string();
    std::string punits = node.child("pressures").attribute("units").as_string();

    // Get the names and stoichiometries of the species that define the reaction
    std::string equation = node.child("equation").text().get();

    // Check if element `temperatures` was provided, if not set default to 25 celsius
    if(temperatures.empty()) temperatures.push_back(25.0);

    // Check if element `pressures` was provided, if not set default to 1 bar
    if(pressures.empty()) pressures.push_back(1.0);

    // Check if temperature units was provided, if not set default to celsius
    if(tunits.empty()) tunits = "celsius";

    // Check if pressure units was provided, if not set default to bar
    if(punits.empty()) punits = "bar";

    // Convert temperatures and pressures to standard units (kelvin and pascal respectively)
    for(auto& x : temperatures) x = units::convert(x, tunits, "kelvin");
    for(auto& x : pressures)    x = units::convert(x, punits, "pascal");

    // Define a lambda function to generate a bilinear interpolator from a vector
    auto bilinear_interpolator = [&](const std::vector<double>& data) -> BilinearInterpolator
    {
        if(data.empty())  return BilinearInterpolator();
        return BilinearInterpolator(temperatures, pressures, data);
    };

    // Define a lambda function to generate a bilinear interpolator of the Gibbs energy of a reaction from its lnk
    auto gibbs_energy_from_lnk = [](const BilinearInterpolator& lnk)
    {
        const double R = universalGasConstant;
        auto f = [=](double T, double P) { return -R*T*lnk(T, P); };
        return BilinearInterpolator(lnk.xCoodinates(), lnk.yCoodinates(), f);
    };

    // Initialise the properties thermodynamic properties of the reaction
    ReactionThermoProperties data;
    data.equation         = equation;
    data.lnk              = bilinear_interpolator(lnk);
    data.gibbs_energy     = gibbs_energy.empty() ? gibbs_energy_from_lnk(data.lnk) : bilinear_interpolator(gibbs_energy);
    data.helmholtz_energy = bilinear_interpolator(helmholtz_energy);
    data.internal_energy  = bilinear_interpolator(internal_energy);
    data.enthalpy         = bilinear_interpolator(enthalpy);
    data.entropy          = bilinear_interpolator(entropy);
    data.volume           = bilinear_interpolator(volume);
    data.heat_capacity    = bilinear_interpolator(heat_capacity);

    return data;
}

auto parseSpeciesThermoProperties(const xml_node& node) -> SpeciesThermoProperties
{
    // Get the data values of the children nodes
    std::vector<double> temperatures     = tofloats(node.child("temperatures").text().get());
    std::vector<double> pressures        = tofloats(node.child("pressures").text().get());
    std::vector<double> gibbs_energy     = tofloats(node.child("G").text().get());
    std::vector<double> helmholtz_energy = tofloats(node.child("A").text().get());
    std::vector<double> internal_energy  = tofloats(node.child("U").text().get());
    std::vector<double> enthalpy         = tofloats(node.child("H").text().get());
    std::vector<double> entropy          = tofloats(node.child("S").text().get());
    std::vector<double> volume           = tofloats(node.child("V").text().get());
    std::vector<double> heat_capacity    = tofloats(node.child("C").text().get());

    // Get the temperature and pressure units
    std::string tunits = node.child("temperatures").attribute("units").as_string();
    std::string punits = node.child("pressures").attribute("units").as_string();

    // Check if element `temperatures` was provided, if not set default to 25 celsius
    if(temperatures.empty()) temperatures.push_back(25.0);

    // Check if element `pressures` was provided, if not set default to 1 bar
    if(pressures.empty()) pressures.push_back(1.0);

    // Check if temperature units was provided, if not set default to celsius
    if(tunits.empty()) tunits = "celsius";

    // Check if pressure units was provided, if not set default to bar
    if(punits.empty()) punits = "bar";

    // Convert temperatures and pressures to standard units
    for(auto& x : temperatures) x = units::convert(x, tunits, "kelvin");
    for(auto& x : pressures)    x = units::convert(x, punits, "pascal");

    // Define a lambda function to generate a bilinear interpolator from a vector
    auto bilinear_interpolator = [&](const std::vector<double>& data) -> BilinearInterpolator
    {
        if(data.empty())  return BilinearInterpolator();
        return BilinearInterpolator(temperatures, pressures, data);
    };

    // Initialise the properties thermodynamic properties of the species
    SpeciesThermoProperties data;
    data.gibbs_energy     = bilinear_interpolator(gibbs_energy);
    data.helmholtz_energy = bilinear_interpolator(helmholtz_energy);
    data.internal_energy  = bilinear_interpolator(internal_energy);
    data.enthalpy         = bilinear_interpolator(enthalpy);
    data.entropy          = bilinear_interpolator(entropy);
    data.volume           = bilinear_interpolator(volume);
    data.heat_capacity    = bilinear_interpolator(heat_capacity);

    return data;
}

auto parseAqueousSpeciesThermoParamsHKF(const xml_node& node) -> Optional<AqueousSpeciesThermoParamsHKF>
{
    const bool emptyGf = node.child("Gf").text().empty();
    const bool emptyHf = node.child("Hf").text().empty();

    if(emptyGf or emptyHf)
        return Optional<AqueousSpeciesThermoParamsHKF>();

    AqueousSpeciesThermoParamsHKF hkf;

    hkf.Gf = node.child("Gf").text().as_double();
    hkf.Hf = node.child("Hf").text().as_double();
    hkf.Sr = node.child("Sr").text().as_double();
    hkf.a1 = node.child("a1").text().as_double();
    hkf.a2 = node.child("a2").text().as_double();
    hkf.a3 = node.child("a3").text().as_double();
    hkf.a4 = node.child("a4").text().as_double();
    hkf.c1 = node.child("c1").text().as_double();
    hkf.c2 = node.child("c2").text().as_double();
    hkf.wref = node.child("wref").text().as_double();

    return hkf;
}

auto parseGaseousSpeciesThermoParamsHKF(const xml_node& node) -> Optional<GaseousSpeciesThermoParamsHKF>
{
    const bool emptyGf = node.child("Gf").text().empty();
    const bool emptyHf = node.child("Hf").text().empty();

    if(emptyGf or emptyHf)
        return Optional<GaseousSpeciesThermoParamsHKF>();

    GaseousSpeciesThermoParamsHKF hkf;

    hkf.Gf = node.child("Gf").text().as_double();
    hkf.Hf = node.child("Hf").text().as_double();
    hkf.Sr = node.child("Sr").text().as_double();
    hkf.a = node.child("a").text().as_double();
    hkf.b = node.child("b").text().as_double();
    hkf.c = node.child("c").text().as_double();
    hkf.Tmax = node.child("Tmax").text().as_double();

    return hkf;
}

auto parseMineralSpeciesThermoParamsHKF(const xml_node& node) -> Optional<MineralSpeciesThermoParamsHKF>
{
    const bool emptyGf = node.child("Gf").text().empty();
    const bool emptyHf = node.child("Hf").text().empty();

    if(emptyGf or emptyHf)
        return Optional<MineralSpeciesThermoParamsHKF>();

    MineralSpeciesThermoParamsHKF hkf;

    hkf.Gf = node.child("Gf").text().as_double();
    hkf.Hf = node.child("Hf").text().as_double();
    hkf.Sr = node.child("Sr").text().as_double();
    hkf.Vr = node.child("Vr").text().as_double();
    hkf.nptrans = node.child("nptrans").text().as_int();
    hkf.Tmax = node.child("Tmax").text().as_double();

    if(hkf.nptrans == 0)
    {
        hkf.a.push_back(node.child("a").text().as_double());
        hkf.b.push_back(node.child("b").text().as_double());
        hkf.c.push_back(node.child("c").text().as_double());
    }
    else
    {
        for(int i = 0; i <= hkf.nptrans; ++i)
        {
            std::stringstream str;
            str << "temperature_range" << i;

            auto temperature_range = node.child(str.str().c_str());

            hkf.a.push_back(temperature_range.child("a").text().as_double());
            hkf.b.push_back(temperature_range.child("b").text().as_double());
            hkf.c.push_back(temperature_range.child("c").text().as_double());

            if(i < hkf.nptrans)
            {
                hkf.Ttr.push_back(temperature_range.child("Ttr").text().as_double());

                const double nan = std::numeric_limits<double>::quiet_NaN();

                const bool empty_Htr = temperature_range.child("Htr").text().empty();
                const bool empty_Vtr = temperature_range.child("Vtr").text().empty();
                const bool empty_dPdTtr = temperature_range.child("dPdTtr").text().empty();

                hkf.Htr.push_back(empty_Htr ? nan : temperature_range.child("Htr").text().as_double());
                hkf.Vtr.push_back(empty_Vtr ? nan : temperature_range.child("Vtr").text().as_double());
                hkf.dPdTtr.push_back(empty_dPdTtr ? nan : temperature_range.child("dPdTtr").text().as_double());
            }
        }
    }

    return hkf;
}

auto parseAqueousSpeciesThermoData(const xml_node& node) -> AqueousSpeciesThermoData
{
    AqueousSpeciesThermoData thermo;

    if(not node.child("properties").empty())
        thermo.properties = parseSpeciesThermoProperties(node);

    if(not node.child("reaction").empty())
        thermo.reaction = parseReactionThermoProperties(node.child("reaction"));

    if(not node.child("hkf").empty())
        thermo.hkf = parseAqueousSpeciesThermoParamsHKF(node.child("hkf"));

    return thermo;
}

auto parseGaseousSpeciesThermoData(const xml_node& node) -> GaseousSpeciesThermoData
{
    GaseousSpeciesThermoData thermo;

    if(not node.child("properties").empty())
        thermo.properties = parseSpeciesThermoProperties(node.child("properties"));

    if(not node.child("reaction").empty())
        thermo.reaction = parseReactionThermoProperties(node.child("reaction"));

    if(not node.child("hkf").empty())
        thermo.hkf = parseGaseousSpeciesThermoParamsHKF(node.child("hkf"));

    return thermo;
}

auto parseMineralSpeciesThermoData(const xml_node& node) -> MineralSpeciesThermoData
{
    MineralSpeciesThermoData thermo;

    if(not node.child("properties").empty())
        thermo.properties = parseSpeciesThermoProperties(node.child("properties"));

    if(not node.child("reaction").empty())
        thermo.reaction = parseReactionThermoProperties(node.child("reaction"));

    if(not node.child("hkf").empty())
        thermo.hkf = parseMineralSpeciesThermoParamsHKF(node.child("hkf"));

    return thermo;
}

template<typename SpeciesMap, typename SpeciesFunction>
auto collectSpecies(const SpeciesMap& map, const SpeciesFunction& fn) -> std::vector<std::string>
{
    std::set<std::string> species;
    for(const auto& entry : map)
        if(fn(entry.second))
            species.insert(entry.first);
    return std::vector<std::string>(species.begin(), species.end());
}

template<typename SpeciesMap>
auto speciesWithElements(const std::vector<std::string>& elements, const SpeciesMap& map) -> std::vector<std::string>
{
    auto f = [&](const Species& species)
    {
        for(std::string element : elements)
            if(species.elementAtoms(element))
                return true;
        return false;
    };

    return collectSpecies(map, f);
}

} // namespace

struct Database::Impl
{
    /// The set of all elements in the database
    ElementMap element_map;

    /// The set of all aqueous species in the database
    AqueousSpeciesMap aqueous_species_map;

    /// The set of all gaseous species in the database
    GaseousSpeciesMap gaseous_species_map;

    /// The set of all mineral species in the database
    MineralSpeciesMap mineral_species_map;

    Impl()
    {}

    Impl(std::string filename)
    {
        parse(filename);
    }

    template<typename Key, typename Value>
    auto collectValues(const std::map<Key, Value>& map) -> std::vector<Value>
    {
        std::vector<Value> species;
        species.reserve(map.size());
        for(const auto& pair : map)
            species.push_back(pair.second);
        return species;
    }

    auto elements() -> std::vector<Element>
    {
        return collectValues(element_map);
    }

    auto aqueousSpecies() -> std::vector<AqueousSpecies>
    {
        return collectValues(aqueous_species_map);
    }

    auto aqueousSpecies(std::string name) const -> const AqueousSpecies&
    {
        if(aqueous_species_map.count(name) == 0)
            errorNonExistentSpecies("aqueous", name);

        return aqueous_species_map.find(name)->second;
    }

    auto gaseousSpecies() -> std::vector<GaseousSpecies>
    {
        return collectValues(gaseous_species_map);
    }

    auto gaseousSpecies(std::string name) const -> const GaseousSpecies&
    {
        if(gaseous_species_map.count(name) == 0)
            errorNonExistentSpecies("gaseous", name);

        return gaseous_species_map.find(name)->second;
    }

    auto mineralSpecies() -> std::vector<MineralSpecies>
    {
        return collectValues(mineral_species_map);
    }

    auto mineralSpecies(std::string name) const -> const MineralSpecies&
    {
        if(mineral_species_map.count(name) == 0)
            errorNonExistentSpecies("mineral", name);

        return mineral_species_map.find(name)->second;
    }

    auto containsAqueousSpecies(std::string species) const -> bool
    {
        return aqueous_species_map.count(species) != 0;
    }

    auto containsGaseousSpecies(std::string species) const -> bool
    {
        return gaseous_species_map.count(species) != 0;
    }

    auto containsMineralSpecies(std::string species) const -> bool
    {
        return mineral_species_map.count(species) != 0;
    }

    auto aqueousSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<std::string>
    {
        return speciesWithElements(elements, aqueous_species_map);
    }

    auto gaseousSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<std::string>
    {
        return speciesWithElements(elements, gaseous_species_map);
    }

    auto mineralSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<std::string>
    {
        return speciesWithElements(elements, mineral_species_map);
    }

    auto parse(std::string filename) -> void
    {
        // Create the XML document
        xml_document doc;

        // Load the xml database file
        auto result = doc.load_file(filename.c_str());

        // Check if the file was correctly loaded
        Assert(result, "Cannot open the database file `" + filename + "`.",
            "The file name or its path might not have been correctly specified.");

        // Access the database node of the database file
        xml_node database = doc.child("database");

        // Read all elements in the database
        for(xml_node node : database.children("element"))
        {
            Element element = parseElement(node);
            element_map[element.name()] = element;
        }

        // Read all species in the database
        for(xml_node node : database.children("species"))
        {
            std::string type = node.child("type").text().get();
            std::string name = node.child("name").text().get();

            if(type == "Aqueous")
            {
                AqueousSpecies species = parseAqueousSpecies(node);
                aqueous_species_map[species.name()] = species;
            }
            else if(type == "Gaseous")
            {
                GaseousSpecies species = parseGaseousSpecies(node);
                gaseous_species_map[species.name()] = species;
            }

            else if(type == "Mineral")
            {
                MineralSpecies species = parseMineralSpecies(node);
                mineral_species_map[species.name()] = species;
            }
            else RuntimeError("Cannot parse the species `" +
                name + "` with type `" + type + "` in the database `" +
                filename + "`.", "The type of the species in unknown.");
        }
    }

    auto parseElement(const xml_node& node) -> Element
    {
        Element element;
        element.setName(node.child("name").text().get());
        element.setMolarMass(node.child("molar_mass").text().as_double());
        return element;
    }

    auto parseElementalFormula(std::string formula) -> std::map<Element, double>
    {
        std::map<Element, double> elements;
        auto words = split(formula, "()");
        for(unsigned i = 0; i < words.size(); i += 2)
        {
            Assert(element_map.count(words[i]),
                "Cannot parse the elemental formula `" + formula + "`.",
                "The element `" + words[i] + "` is not in the database.");
            elements.emplace(element_map.at(words[i]), tofloat(words[i + 1]));
        }
        return elements;
    }

    auto parseSpecies(const xml_node& node) -> Species
    {
        // The species instance
        Species species;

        // Set the name of the species
        species.setName(node.child("name").text().get());

        // Set the chemical formula of the species
        species.setFormula(node.child("formula").text().get());

        // Set the elements of the species
        species.setElements(parseElementalFormula(node.child("elements").text().get()));

        // Set the molar mass of the species
        if(not node.child("molar_mass").text().empty())
        {
            const auto value = node.child("molar_mass").text().as_double();
            const auto units = node.child("molar_mass").attribute("units").as_string();
            species.setMolarMass(units::convert(value, units, "kg/mol"));
        }

        return species;
    }

    auto parseAqueousSpecies(const xml_node& node) -> AqueousSpecies
    {
        // The aqueous species instance
        AqueousSpecies species = parseSpecies(node);

        // Set the elemental charge of the species
        species.setCharge(node.child("charge").text().as_double());

        // Parse the complex formula of the aqueous species (if any)
        species.setDissociation(parseDissociation(node.child("dissociation").text().get()));

        // Parse the thermodynamic data of the aqueous species
        species.setThermoData(parseAqueousSpeciesThermoData(node.child("thermo")));

        return species;
    }

    auto parseGaseousSpecies(const xml_node& node) -> GaseousSpecies
    {
        // The gaseous species instance
        GaseousSpecies species = parseSpecies(node);

        // Parse the thermodynamic data of the gaseous species
        species.setThermoData(parseGaseousSpeciesThermoData(node.child("thermo")));

        return species;
    }

    auto parseMineralSpecies(const xml_node& node) -> MineralSpecies
    {
        // The mineral species instance
        MineralSpecies species = parseSpecies(node);

        // Parse the thermodynamic data of the mineral species
        species.setThermoData(parseMineralSpeciesThermoData(node.child("thermo")));

        return species;
    }
};

Database::Database()
: pimpl(new Impl())
{}

Database::Database(std::string filename)
: pimpl(new Impl(filename))
{}

auto Database::elements() -> std::vector<Element>
{
    return pimpl->elements();
}

auto Database::aqueousSpecies() -> std::vector<AqueousSpecies>
{
    return pimpl->aqueousSpecies();
}

auto Database::aqueousSpecies(std::string name) const -> const AqueousSpecies&
{
    return pimpl->aqueousSpecies(name);
}

auto Database::gaseousSpecies() -> std::vector<GaseousSpecies>
{
    return pimpl->gaseousSpecies();
}

auto Database::gaseousSpecies(std::string name) const -> const GaseousSpecies&
{
    return pimpl->gaseousSpecies(name);
}

auto Database::mineralSpecies() -> std::vector<MineralSpecies>
{
    return pimpl->mineralSpecies();
}

auto Database::mineralSpecies(std::string name) const -> const MineralSpecies&
{
    return pimpl->mineralSpecies(name);
}

auto Database::containsAqueousSpecies(std::string species) const -> bool
{
    return pimpl->containsAqueousSpecies(species);
}

auto Database::containsGaseousSpecies(std::string species) const -> bool
{
    return pimpl->containsGaseousSpecies(species);
}

auto Database::containsMineralSpecies(std::string species) const -> bool
{
    return pimpl->containsMineralSpecies(species);
}

auto Database::aqueousSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<std::string>
{
    return pimpl->aqueousSpeciesWithElements(elements);
}

auto Database::gaseousSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<std::string>
{
    return pimpl->gaseousSpeciesWithElements(elements);
}

auto Database::mineralSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<std::string>
{
    return pimpl->mineralSpeciesWithElements(elements);
}

} // namespace Reaktor

