// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Thermodynamics/Databases/DatabaseUtils.hpp>
#include <Reaktoro/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/GaseousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/MineralSpecies.hpp>

// miniz includes
#include <miniz/zip_file.hpp>

// pugixml includes
#include <pugixml.hpp>
using namespace pugi;

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

auto parseReactionInterpolatedThermoProperties(const xml_node& node) -> ReactionThermoInterpolatedProperties
{
    // Get the data values of the children nodes
    std::vector<double> temperatures     = tofloats(node.child("Temperatures").text().get());
    std::vector<double> pressures        = tofloats(node.child("Pressures").text().get());
    std::vector<double> pk               = tofloats(node.child("pk").text().get());
    std::vector<double> lnk              = tofloats(node.child("lnk").text().get());
    std::vector<double> logk             = tofloats(node.child("logk").text().get());
    std::vector<double> gibbs_energy     = tofloats(node.child("G").text().get());
    std::vector<double> helmholtz_energy = tofloats(node.child("A").text().get());
    std::vector<double> internal_energy  = tofloats(node.child("U").text().get());
    std::vector<double> enthalpy         = tofloats(node.child("H").text().get());
    std::vector<double> entropy          = tofloats(node.child("S").text().get());
    std::vector<double> volume           = tofloats(node.child("V").text().get());
    std::vector<double> heat_capacity_cp = tofloats(node.child("Cp").text().get());
    std::vector<double> heat_capacity_cv = tofloats(node.child("Cv").text().get());

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
    std::string tunits = node.child("Temperatures").attribute("units").as_string();
    std::string punits = node.child("Pressures").attribute("units").as_string();

    // Get the names and stoichiometries of the species that define the reaction
    std::string equation = node.child("Equation").text().get();

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

    // Initialize the properties thermodynamic properties of the reaction
    ReactionThermoInterpolatedProperties data;
    data.equation         = equation;
    data.lnk              = bilinear_interpolator(lnk);
    data.gibbs_energy     = gibbs_energy.empty() ? gibbs_energy_from_lnk(data.lnk) : bilinear_interpolator(gibbs_energy);
    data.helmholtz_energy = bilinear_interpolator(helmholtz_energy);
    data.internal_energy  = bilinear_interpolator(internal_energy);
    data.enthalpy         = bilinear_interpolator(enthalpy);
    data.entropy          = bilinear_interpolator(entropy);
    data.volume           = bilinear_interpolator(volume);
    data.heat_capacity_cp = bilinear_interpolator(heat_capacity_cp);
    data.heat_capacity_cv = bilinear_interpolator(heat_capacity_cv);

    return data;
}

auto parseSpeciesInterpolatedThermoProperties(const xml_node& node) -> SpeciesThermoInterpolatedProperties
{
    // Get the data values of the children nodes
    std::vector<double> temperatures     = tofloats(node.child("Temperatures").text().get());
    std::vector<double> pressures        = tofloats(node.child("Pressures").text().get());
    std::vector<double> gibbs_energy     = tofloats(node.child("G").text().get());
    std::vector<double> helmholtz_energy = tofloats(node.child("A").text().get());
    std::vector<double> internal_energy  = tofloats(node.child("U").text().get());
    std::vector<double> enthalpy         = tofloats(node.child("H").text().get());
    std::vector<double> entropy          = tofloats(node.child("S").text().get());
    std::vector<double> volume           = tofloats(node.child("V").text().get());
    std::vector<double> heat_capacity_cp = tofloats(node.child("Cp").text().get());
    std::vector<double> heat_capacity_cv = tofloats(node.child("Cv").text().get());

    // Get the temperature and pressure units
    std::string tunits = node.child("Temperatures").attribute("units").as_string();
    std::string punits = node.child("Pressures").attribute("units").as_string();

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

    // Initialize the properties thermodynamic properties of the species
    SpeciesThermoInterpolatedProperties data;
    data.gibbs_energy     = bilinear_interpolator(gibbs_energy);
    data.helmholtz_energy = bilinear_interpolator(helmholtz_energy);
    data.internal_energy  = bilinear_interpolator(internal_energy);
    data.enthalpy         = bilinear_interpolator(enthalpy);
    data.entropy          = bilinear_interpolator(entropy);
    data.volume           = bilinear_interpolator(volume);
    data.heat_capacity_cp = bilinear_interpolator(heat_capacity_cp);
    data.heat_capacity_cv = bilinear_interpolator(heat_capacity_cv);

    return data;
}

auto parseAqueousSpeciesThermoParamsHKF(const xml_node& node) -> Optional<AqueousSpeciesThermoParamsHKF>
{
    const bool emptyGf = node.child("Gf").text().empty();
    const bool emptyHf = node.child("Hf").text().empty();

    if(emptyGf || emptyHf)
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

    if(emptyGf || emptyHf)
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

    if(emptyGf || emptyHf)
        return Optional<MineralSpeciesThermoParamsHKF>();

    MineralSpeciesThermoParamsHKF hkf;

    hkf.Gf = node.child("Gf").text().as_double();
    hkf.Hf = node.child("Hf").text().as_double();
    hkf.Sr = node.child("Sr").text().as_double();
    hkf.Vr = node.child("Vr").text().as_double();
    hkf.nptrans = node.child("NumPhaseTrans").text().as_int();
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
            str << "TemperatureRange" << i;

            auto temperature_range = node.child(str.str().c_str());

            hkf.a.push_back(temperature_range.child("a").text().as_double());
            hkf.b.push_back(temperature_range.child("b").text().as_double());
            hkf.c.push_back(temperature_range.child("c").text().as_double());

            if(i < hkf.nptrans)
            {
                hkf.Ttr.push_back(temperature_range.child("Ttr").text().as_double());

                // Check if deltaH, deltaV and dPdT for phase transition is available
                const bool empty_Htr = temperature_range.child("Htr").text().empty();
                const bool empty_Vtr = temperature_range.child("Vtr").text().empty();
                const bool empty_dPdTtr = temperature_range.child("dPdTtr").text().empty();

                // Set zero the non-available transition values
                hkf.Htr.push_back(empty_Htr ? 0.0 : temperature_range.child("Htr").text().as_double());
                hkf.Vtr.push_back(empty_Vtr ? 0.0 : temperature_range.child("Vtr").text().as_double());
                hkf.dPdTtr.push_back(empty_dPdTtr ? 0.0 : temperature_range.child("dPdTtr").text().as_double());
            }
        }
    }

    return hkf;
}

auto parseAqueousSpeciesThermoData(const xml_node& node) -> AqueousSpeciesThermoData
{
    AqueousSpeciesThermoData thermo;

    if(!node.child("Properties").empty())
        thermo.properties = parseSpeciesInterpolatedThermoProperties(node.child("Properties"));

    if(!node.child("Reaction").empty())
        thermo.reaction = parseReactionInterpolatedThermoProperties(node.child("Reaction"));

    if(!node.child("HKF").empty())
        thermo.hkf = parseAqueousSpeciesThermoParamsHKF(node.child("HKF"));

    return thermo;
}

auto parseGaseousSpeciesThermoData(const xml_node& node) -> GaseousSpeciesThermoData
{
    GaseousSpeciesThermoData thermo;

    if(!node.child("Properties").empty())
        thermo.properties = parseSpeciesInterpolatedThermoProperties(node.child("Properties"));

    if(!node.child("Reaction").empty())
        thermo.reaction = parseReactionInterpolatedThermoProperties(node.child("Reaction"));

    if(!node.child("HKF").empty())
        thermo.hkf = parseGaseousSpeciesThermoParamsHKF(node.child("HKF"));

    return thermo;
}

auto parseMineralSpeciesThermoData(const xml_node& node) -> MineralSpeciesThermoData
{
    MineralSpeciesThermoData thermo;

    if(!node.child("Properties").empty())
        thermo.properties = parseSpeciesInterpolatedThermoProperties(node.child("Properties"));

    if(!node.child("Reaction").empty())
        thermo.reaction = parseReactionInterpolatedThermoProperties(node.child("Reaction"));

    if(!node.child("HKF").empty())
        thermo.hkf = parseMineralSpeciesThermoParamsHKF(node.child("HKF"));

    return thermo;
}

template<typename SpeciesType, typename SpeciesFunction>
auto collectSpecies(const std::map<std::string, SpeciesType>& map, const SpeciesFunction& fn) -> std::vector<SpeciesType>
{
    std::set<SpeciesType> species;
    for(const auto& entry : map)
        if(fn(entry.second))
            species.insert(entry.second);
    return std::vector<SpeciesType>(species.begin(), species.end());
}

template<typename SpeciesType>
auto speciesWithElements(const std::vector<std::string>& elements, const std::map<std::string, SpeciesType>& map) -> std::vector<SpeciesType>
{
    auto f = [&](const SpeciesType& species)
    {
        // Check if the given species has all chemical elements in the list of elements (ignore charge Z)
        for(auto pair : species.elements())
            if(pair.first.name() != "Z" && !contained(pair.first.name(), elements))
                return false;
        return true;
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
        // Create the XML document
        xml_document doc;

        // Load the xml database file
        auto result = doc.load_file(filename.c_str());

        // Check if result is not ok, and then try a built-in database with same name
        if(!result)
        {
            // Search for a built-in database
            std::string builtin = database(filename);

            // If not empty, use the built-in database to create the xml doc
            if(!builtin.empty()) result = doc.load(builtin.c_str());
        }

        // Ensure either a database file path was correctly given, or a built-in database
        if(!result)
        {
            std::string names;
            for(auto const& s : databases()) names += s + " ";
            RuntimeError("Could not initialize the Database instance with given database name `" + filename + "`.",
                "This name either points to a non-existent database file, or it is not one of the "
                "built-in database files in Reaktoro. The built-in databases are: " + names + ".");
        }

        // Parse the xml document
        parse(doc, filename);
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

    auto addElement(const Element& element) -> void
    {
    	element_map.insert({element.name(), element});
    }

    auto addAqueousSpecies(const AqueousSpecies& species) -> void
    {
    	aqueous_species_map.insert({species.name(), species});
    }

    auto addGaseousSpecies(const GaseousSpecies& species) -> void
    {
    	gaseous_species_map.insert({species.name(), species});
    }

    auto addMineralSpecies(const MineralSpecies& species) -> void
    {
    	mineral_species_map.insert({species.name(), species});
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

    auto aqueousSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<AqueousSpecies>
    {
        return speciesWithElements(elements, aqueous_species_map);
    }

    auto gaseousSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<GaseousSpecies>
    {
        return speciesWithElements(elements, gaseous_species_map);
    }

    auto mineralSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<MineralSpecies>
    {
        return speciesWithElements(elements, mineral_species_map);
    }

    auto parse(const xml_document& doc, std::string databasename) -> void
    {
        // Access the database node of the database file
        xml_node database = doc.child("Database");

        // Read all elements in the database
        for(xml_node node : database.children("Element"))
        {
            Element element = parseElement(node);
            element_map[element.name()] = element;
        }

        // Read all species in the database
        for(xml_node node : database.children("Species"))
        {
            std::string type = node.child("Type").text().get();
            std::string name = node.child("Name").text().get();

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
            else RuntimeError("Could not parse the species `" +
                name + "` with type `" + type + "` in the database `" +
                databasename + "`.", "The type of the species in unknown.");
        }
    }

    auto parseElement(const xml_node& node) -> Element
    {
        Element element;
        element.setName(node.child("Name").text().get());
        element.setMolarMass(node.child("MolarMass").text().as_double() * 1e-3); // convert from g/mol to kg/mol
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
        species.setName(node.child("Name").text().get());

        // Set the chemical formula of the species
        species.setFormula(node.child("Formula").text().get());

        // Set the elements of the species
        species.setElements(parseElementalFormula(node.child("Elements").text().get()));

        return species;
    }

    auto parseAqueousSpecies(const xml_node& node) -> AqueousSpecies
    {
        // The aqueous species instance
        AqueousSpecies species = parseSpecies(node);

        // Set the elemental charge of the species
        species.setCharge(node.child("Charge").text().as_double());

        // Parse the complex formula of the aqueous species (if any)
        species.setDissociation(parseDissociation(node.child("Dissociation").text().get()));

        // Parse the thermodynamic data of the aqueous species
        species.setThermoData(parseAqueousSpeciesThermoData(node.child("Thermo")));

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

    auto parseGaseousSpecies(const xml_node& node) -> GaseousSpecies
    {
        // The gaseous species instance
        GaseousSpecies species = parseSpecies(node);

        // Set the critical temperature of the gaseous species (in units of K)
        if(!node.child("CriticalTemperature").empty())
            species.setCriticalTemperature(node.child("CriticalTemperature").text().as_double());

        // Set the critical pressure of the gaseous species (in units of Pa)
        if(!node.child("CriticalPressure").empty())
            species.setCriticalPressure(node.child("CriticalPressure").text().as_double() * 1e5); // convert from bar to Pa

        // Set the acentric factor of the gaseous species
        if(!node.child("AcentricFactor").empty())
            species.setAcentricFactor(node.child("AcentricFactor").text().as_double());

        // Parse the thermodynamic data of the gaseous species
        species.setThermoData(parseGaseousSpeciesThermoData(node.child("Thermo")));

        return species;
    }

    auto parseMineralSpecies(const xml_node& node) -> MineralSpecies
    {
        // The mineral species instance
        MineralSpecies species = parseSpecies(node);

        // Parse the thermodynamic data of the mineral species
        species.setThermoData(parseMineralSpeciesThermoData(node.child("Thermo")));

        return species;
    }
};

Database::Database()
: pimpl(new Impl())
{}

Database::Database(std::string filename)
: pimpl(new Impl(filename))
{}

auto Database::addElement(const Element& element) -> void
{
    pimpl->addElement(element);
}

auto Database::addAqueousSpecies(const AqueousSpecies& species) -> void
{
    pimpl->addAqueousSpecies(species);
}

auto Database::addGaseousSpecies(const GaseousSpecies& species) -> void
{
    pimpl->addGaseousSpecies(species);
}

auto Database::addMineralSpecies(const MineralSpecies& species) -> void
{
    pimpl->addMineralSpecies(species);
}

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

auto Database::aqueousSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<AqueousSpecies>
{
    return pimpl->aqueousSpeciesWithElements(elements);
}

auto Database::gaseousSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<GaseousSpecies>
{
    return pimpl->gaseousSpeciesWithElements(elements);
}

auto Database::mineralSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<MineralSpecies>
{
    return pimpl->mineralSpeciesWithElements(elements);
}

} // namespace Reaktoro

