// Reaktoro is a C++ library for computational reaction modelling.
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

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Optional.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Thermodynamics/Species/GeneralSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/GaseousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/MineralSpecies.hpp>

// pugixml includes
#include <yaml-cpp/include/yaml-cpp/yaml.h>

namespace Reaktoro {
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

auto extractValue(const YAML::Node& node) -> double
{
    if(node.Type() == YAML::NodeType::Null) return 0.0;

    if(node["value"])
        return node["value"].as<double>();
    else
        return node.as<double>();
}

auto extractValues(const YAML::Node& node) -> std::vector<double>
{
    if(node.Type() == YAML::NodeType::Null) return {};

    if(node["values"])
        return tofloats(node["values"].as<std::string>());
    else
        return tofloats(node.as<std::string>());
}

auto extractUnits(const YAML::Node& node, std::string default_units) -> std::string
{
    if(node.Type() == YAML::NodeType::Null) return default_units;

    if(node["units"])
        return node["units"].as<std::string>();
    else
        return default_units;
}

auto parseReactionThermoProperties(const YAML::Node& node) -> ReactionThermoProperties
{
    // Assert a reaction equation was specified
    Assert(node["Equation"], "Could not initialize the database.",
        "Expecting a reaction equation via the `Equation` keyword in "
        "the `Reaction` data block.")

    // Get the data values of the children nodes
    std::vector<double> temperatures     = extractValues(node["Temperatures"]);
    std::vector<double> pressures        = extractValues(node["Pressures"]);
    std::vector<double> pk               = extractValues(node["pk"]);
    std::vector<double> lnk              = extractValues(node["lnk"]);
    std::vector<double> logk             = extractValues(node["logk"]);
    std::vector<double> gibbs_energy     = extractValues(node["G"]);
    std::vector<double> helmholtz_energy = extractValues(node["A"]);
    std::vector<double> internal_energy  = extractValues(node["U"]);
    std::vector<double> enthalpy         = extractValues(node["H"]);
    std::vector<double> entropy          = extractValues(node["S"]);
    std::vector<double> volume           = extractValues(node["V"]);
    std::vector<double> heat_capacity    = extractValues(node["C"]);

    // Convert `pk` to `lnk`, where `pk = -log(k) = -ln(k)/ln(10)`
    const double ln_10 = std::log(10.0);
    if(!pk.empty() && lnk.empty())
    {
        lnk.resize(pk.size());
        for(unsigned i = 0; i < pk.size(); ++i)
            lnk[i] = -pk[i] * ln_10;
    }

    // Convert `logk` to `lnk`, where `log(k) = ln(k)/ln(10)`
    if(!logk.empty() && lnk.empty())
    {
        lnk.resize(logk.size());
        for(unsigned i = 0; i < logk.size(); ++i)
            lnk[i] = logk[i] * ln_10;
    }

    // Get the temperature and pressure units
    std::string tunits = extractUnits(node["Temperatures"], "celsius");
    std::string punits = extractUnits(node["Pressures"], "bar");

    // Get the names and stoichiometries of the species that define the reaction
    std::string equation = node["Equation"].as<std::string>();

    // Check if element `temperatures` was provided, if not set default to 25 celsius
    if(temperatures.empty()) temperatures.push_back(25.0);

    // Check if element `pressures` was provided, if not set default to 1 bar
    if(pressures.empty()) pressures.push_back(1.0);

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

    // Initialize the properties thermodynamic properties of the reaction
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

auto parseSpeciesThermoProperties(const YAML::Node& node) -> SpeciesThermoProperties
{
    // Get the data values of the children nodes
    std::vector<double> temperatures     = extractValues(node["Temperatures"]);
    std::vector<double> pressures        = extractValues(node["Pressures"]);
    std::vector<double> gibbs_energy     = extractValues(node["G"]);
    std::vector<double> helmholtz_energy = extractValues(node["A"]);
    std::vector<double> internal_energy  = extractValues(node["U"]);
    std::vector<double> enthalpy         = extractValues(node["H"]);
    std::vector<double> entropy          = extractValues(node["S"]);
    std::vector<double> volume           = extractValues(node["V"]);
    std::vector<double> heat_capacity    = extractValues(node["C"]);

    // Get the temperature and pressure units
    std::string tunits = extractUnits(node["Temperatures"], "celsius");
    std::string punits = extractUnits(node["Pressures"], "bar");

    // Check if element `temperatures` was provided, if not set default to 25 celsius
    if(temperatures.empty()) temperatures.push_back(25.0);

    // Check if element `pressures` was provided, if not set default to 1 bar
    if(pressures.empty()) pressures.push_back(1.0);

    // Convert temperatures and pressures to standard units
    for(auto& x : temperatures) x = units::convert(x, tunits, "kelvin");
    for(auto& x : pressures)    x = units::convert(x, punits, "pascal");

    // Define a lambda function to generate a bilinear interpolator from a vector
    auto bilinear_interpolator = [&](const std::vector<double>& data) -> BilinearInterpolator
    {
        if(data.empty())  return BilinearInterpolator();
        return BilinearInterpolator(temperatures, pressures, data);
    };

    // Initialize the properties thermodynamic properties of the species
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

auto parseAqueousSpeciesThermoParamsHKF(const YAML::Node& node) -> Optional<AqueousSpeciesThermoParamsHKF>
{
    if(!node["Gf"] || !node["Hf"])
        return Optional<AqueousSpeciesThermoParamsHKF>();

    AqueousSpeciesThermoParamsHKF hkf;

    hkf.Gf = extractValue(node["Gf"]);
    hkf.Hf = extractValue(node["Hf"]);
    hkf.Sr = extractValue(node["Sr"]);
    hkf.a1 = extractValue(node["a1"]);
    hkf.a2 = extractValue(node["a2"]);
    hkf.a3 = extractValue(node["a3"]);
    hkf.a4 = extractValue(node["a4"]);
    hkf.c1 = extractValue(node["c1"]);
    hkf.c2 = extractValue(node["c2"]);
    hkf.wref = extractValue(node["wref"]);

    return hkf;
}

auto parseGaseousSpeciesThermoParamsHKF(const YAML::Node& node) -> Optional<GaseousSpeciesThermoParamsHKF>
{
    if(!node["Gf"] || !node["Hf"])
        return Optional<GaseousSpeciesThermoParamsHKF>();

    GaseousSpeciesThermoParamsHKF hkf;

    hkf.Gf = extractValue(node["Gf"]);
    hkf.Hf = extractValue(node["Hf"]);
    hkf.Sr = extractValue(node["Sr"]);
    hkf.a = extractValue(node["a"]);
    hkf.b = extractValue(node["b"]);
    hkf.c = extractValue(node["c"]);
    hkf.Tmax = extractValue(node["Tmax"]);

    return hkf;
}

auto parseMineralSpeciesThermoParamsHKF(const YAML::Node& node) -> Optional<MineralSpeciesThermoParamsHKF>
{
    if(!node["Gf"] || !node["Hf"])
        return Optional<MineralSpeciesThermoParamsHKF>();

    MineralSpeciesThermoParamsHKF hkf;

    hkf.Gf = extractValue(node["Gf"]);
    hkf.Hf = extractValue(node["Hf"]);
    hkf.Sr = extractValue(node["Sr"]);
    hkf.Vr = extractValue(node["Vr"]);
    hkf.nptrans = node["NumPhaseTrans"].as<int>();
    hkf.Tmax = extractValue(node["Tmax"]);

    if(hkf.nptrans == 0)
    {
        hkf.a.push_back(extractValue(node["a"]));
        hkf.b.push_back(extractValue(node["b"]));
        hkf.c.push_back(extractValue(node["c"]));
    }
    else
    {
        for(int i = 0; i <= hkf.nptrans; ++i)
        {
            auto tag = "TemperatureRange" + std::to_string(i);

            auto temperature_range = node[tag];

            hkf.a.push_back(extractValue(temperature_range["a"]));
            hkf.b.push_back(extractValue(temperature_range["b"]));
            hkf.c.push_back(extractValue(temperature_range["c"]));

            if(i < hkf.nptrans)
            {
                hkf.Ttr.push_back(extractValue(temperature_range["Ttr"]));

                const double nan = std::numeric_limits<double>::quiet_NaN();

                const bool empty_Htr = !temperature_range["Htr"];
                const bool empty_Vtr = !temperature_range["Vtr"];
                const bool empty_dPdTtr = !temperature_range["dPdTtr"];

                hkf.Htr.push_back(empty_Htr ? nan : extractValue(temperature_range["Htr"]));
                hkf.Vtr.push_back(empty_Vtr ? nan : extractValue(temperature_range["Vtr"]));
                hkf.dPdTtr.push_back(empty_dPdTtr ? nan : extractValue(temperature_range["dPdTtr"]));
            }
        }
    }

    return hkf;
}

auto parseAqueousSpeciesThermoData(const YAML::Node& node) -> AqueousSpeciesThermoData
{
    AqueousSpeciesThermoData thermo;

    if(node["Properties"])
        thermo.properties = parseSpeciesThermoProperties(node["Properties"]);

    if(node["Reaction"])
        thermo.reaction = parseReactionThermoProperties(node["Reaction"]);

    if(node["HKF"])
        thermo.hkf = parseAqueousSpeciesThermoParamsHKF(node["HKF"]);

    return thermo;
}

auto parseGaseousSpeciesThermoData(const YAML::Node& node) -> GaseousSpeciesThermoData
{
    GaseousSpeciesThermoData thermo;

    if(node["Properties"])
        thermo.properties = parseSpeciesThermoProperties(node["Properties"]);

    if(node["Reaction"])
        thermo.reaction = parseReactionThermoProperties(node["Reaction"]);

    if(node["HKF"])
        thermo.hkf = parseGaseousSpeciesThermoParamsHKF(node["HKF"]);

    return thermo;
}

auto parseMineralSpeciesThermoData(const YAML::Node& node) -> MineralSpeciesThermoData
{
    MineralSpeciesThermoData thermo;

    if(node["Properties"])
        thermo.properties = parseSpeciesThermoProperties(node["Properties"]);

    if(node["Reaction"])
        thermo.reaction = parseReactionThermoProperties(node["Reaction"]);

    if(node["HKF"])
        thermo.hkf = parseMineralSpeciesThermoParamsHKF(node["HKF"]);

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
    auto f = [&](const GeneralSpecies& species)
    {
        for(std::string element : elements)
            if(species.elementCoefficient(element))
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
        // Create the YAML document
        YAML::Node doc = YAML::LoadFile(filename);

        // Check if the file was correctly loaded
        Assert(doc, "Could not initialize the database.",
            "The given file name `" + filename + "` could not be found, or parsed. "
            "Ensure the relative path to the file, with respect to the application, is also specified.");

        // Read all elements in the database
        for(const YAML::Node& node : doc["Elements"])
        {
            Element element = parseElement(node);
            element_map[element.name()] = element;
        }

        // Read all species in the database
        for(const YAML::Node node : doc["Species"])
        {
            std::string type = node["Type"].as<std::string>();
            std::string name = node["Name"].as<std::string>();

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
            else RuntimeError("Could not initialize the database.",
                "Could not parse the species `" + name + "` with unknown type `" +
                type + "` in the database `" + filename + "`.");
        }
    }

    auto parseElement(const YAML::Node& node) -> Element
    {
        Element element;
        element.setName(node["Name"].as<std::string>());
        element.setMolarMass(extractValue(node["MolarMass"]));
        return element;
    }

    auto parseElementalFormula(std::string formula) -> std::map<Element, double>
    {
        std::map<Element, double> elements;
        auto words = split(formula, "()");
        for(unsigned i = 0; i < words.size(); i += 2)
        {
            Assert(element_map.count(words[i]),
                "Could not initialize the database.",
                "Could not parse the elemental formula `" + formula + "` because "
                "the element `" + words[i] + "` is not in the database.");
            elements.emplace(element_map.at(words[i]), tofloat(words[i + 1]));
        }
        return elements;
    }

    auto parseSpecies(const YAML::Node& node) -> GeneralSpecies
    {
        // The species instance
        GeneralSpecies species;

        // Set the name of the species
        species.setName(node["Name"].as<std::string>());

        // Set the chemical formula of the species
        species.setFormula(node["Formula"].as<std::string>());

        // Set the elements of the species
        species.setElements(parseElementalFormula(node["Elements"].as<std::string>()));

        // Set the molar mass of the species
        if(node["MolarMass"])
        {
            const auto value = extractValue(node["MolarMass"]);
            const auto units = extractUnits(node["MolarMass"], "g/mol");
            species.setMolarMass(units::convert(value, units, "kg/mol"));
        }

        return species;
    }

    auto parseAqueousSpecies(const YAML::Node& node) -> AqueousSpecies
    {
        // The aqueous species instance
        AqueousSpecies species = parseSpecies(node);

        // Set the elemental charge of the species
        species.setCharge(node["Charge"].as<double>());

        // Parse the complex formula of the aqueous species (if any)
        if(node["Dissociation"])
            species.setDissociation(parseDissociation(node["Dissociation"].as<std::string>()));

        // Parse the thermodynamic data of the aqueous species
        if(node["Thermo"])
            species.setThermoData(parseAqueousSpeciesThermoData(node["Thermo"]));

        // Update the list of elements of the aqueous species if it is electrically charged
        if(species.charge())
        {
            Element charge;
            charge.setName("Z");
            auto elements = species.elements();
            elements.emplace(charge, species.charge());
            species.setElements(elements);
        }

        return species;
    }

    auto parseGaseousSpecies(const YAML::Node& node) -> GaseousSpecies
    {
        // The gaseous species instance
        GaseousSpecies species = parseSpecies(node);

        // Parse the thermodynamic data of the gaseous species
        species.setThermoData(parseGaseousSpeciesThermoData(node["Thermo"]));

        return species;
    }

    auto parseMineralSpecies(const YAML::Node& node) -> MineralSpecies
    {
        // The mineral species instance
        MineralSpecies species = parseSpecies(node);

        // Parse the thermodynamic data of the mineral species
        species.setThermoData(parseMineralSpeciesThermoData(node["Thermo"]));

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

} // namespace Reaktoro

