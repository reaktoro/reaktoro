/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "Database.hpp"

// C++ includes
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Common/Constants.hpp>
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Common/Units.hpp>
#include <Reaktor/Species/AqueousSpecies.hpp>
#include <Reaktor/Species/GaseousSpecies.hpp>
#include <Reaktor/Species/GeneralSpecies.hpp>
#include <Reaktor/Species/MineralSpecies.hpp>
#include <Reaktor/Utils/SetUtils.hpp>
#include <Reaktor/Utils/StringUtils.hpp>

// pugixml includes
#include <Reaktor/External/pugixml/pugixml.hpp>
using namespace pugi;

namespace Reaktor {
namespace internal {

auto nonexistentSpeciesError(const std::string& type, const std::string& name) -> void
{
    Exception exception;
    exception.error << "Cannot get an instance of the " << type << " species named " << name << " in the database.";
    exception.reason << "There is no such species in the database.";
    raise(exception);
}

BilinearInterpolator
deltaGibssFromLnK(const BilinearInterpolator& lnK)
{
    const double R = universalGasConstant;

    auto f = [=](double T, double P)
    {
        return -R*T*lnK(T, P);
    };

    return BilinearInterpolator(lnK.xCoodinates(), lnK.yCoodinates(), f);
}

BilinearInterpolator
deltaGibssFromLogK(const BilinearInterpolator& logK)
{
    const double R = universalGasConstant;
    const double ln10 = std::log(10.0);

    auto f = [=](double T, double P)
    {
        return -R*T*logK(T, P) * ln10;
    };

    return BilinearInterpolator(logK.xCoodinates(), logK.yCoodinates(), f);
}

BilinearInterpolator
deltaGibssFromPK(const BilinearInterpolator& pK)
{
    const double R = universalGasConstant;
    const double ln10 = std::log(10.0);

    auto f = [=](double T, double P)
    {
        return R*T*pK(T, P) * ln10;
    };

    return BilinearInterpolator(pK.xCoodinates(), pK.yCoodinates(), f);
}

auto parseThermoDataInterpolated(const xml_node& node) -> ThermoDataSpecies
{
    // Get the data values of the children nodes
    std::vector<double> temperatures = tofloats(node.child("temperatures").text().get());
    std::vector<double> pressures    = tofloats(node.child("pressures").text().get());
    std::vector<double> gibbs        = tofloats(node.child("gibbs").text().get());
    std::vector<double> enthalpy     = tofloats(node.child("enthalpy").text().get());
    std::vector<double> entropy      = tofloats(node.child("entropy").text().get());
    std::vector<double> volume       = tofloats(node.child("volume").text().get());
    std::vector<double> cp           = tofloats(node.child("cp").text().get());

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

    auto bilinear_interpolator = [&](const std::vector<double>& data)
    {
        return data.empty() ? BilinearInterpolator() : BilinearInterpolator(temperatures, pressures, data);
    };

    // Initialise the interpolated thermodynamic properties of the species
    ThermoDataSpecies data;

    data.gibbs    = bilinear_interpolator(gibbs);
    data.enthalpy = bilinear_interpolator(enthalpy);
    data.entropy  = bilinear_interpolator(entropy);
    data.volume   = bilinear_interpolator(volume);
    data.cp       = bilinear_interpolator(cp);

    return data;
}

auto parseThermoDataReaction(const xml_node& node) -> ThermoDataReaction
{
    // Get the data values of the children nodes
    std::vector<double> temperatures = tofloats(node.child("temperatures").text().get());
    std::vector<double> pressures    = tofloats(node.child("pressures").text().get());
    std::vector<double> pK           = tofloats(node.child("pK").text().get());
    std::vector<double> lnK          = tofloats(node.child("lnK").text().get());
    std::vector<double> logK         = tofloats(node.child("logK").text().get());
    std::vector<double> gibbs        = tofloats(node.child("gibbs").text().get());
    std::vector<double> enthalpy     = tofloats(node.child("enthalpy").text().get());
    std::vector<double> entropy      = tofloats(node.child("entropy").text().get());
    std::vector<double> volume       = tofloats(node.child("volume").text().get());

    // Get the temperature and pressure units
    std::string tunits = node.child("temperatures").attribute("units").as_string();
    std::string punits = node.child("pressures").attribute("units").as_string();

    // Get the names and stoichiometries of the species that define the reaction
    std::vector<std::string> species = split(node.child("equation").child("species").text().get());
    std::vector<double> stoichiometries = tofloats(node.child("equation").child("stoichiometries").text().get());

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

    auto bilinear_interpolator = [&](const std::vector<double>& data)
    {
        return data.empty() ? BilinearInterpolator() : BilinearInterpolator(temperatures, pressures, data);
    };

    // Initialise the interpolated thermodynamic properties of the reaction
    ThermoDataReaction data;

    // Set the reaction equation for which the standard molar thermodynamic properties are provided
    data.equation = ReactionEquation(species, stoichiometries);

    // Set the equilibrium constants and thermodynamic properties of the reaction
    data.pK       = bilinear_interpolator(pK);
    data.lnK      = bilinear_interpolator(lnK);
    data.logK     = bilinear_interpolator(logK);
    data.gibbs    = bilinear_interpolator(gibbs);
    data.enthalpy = bilinear_interpolator(enthalpy);
    data.entropy  = bilinear_interpolator(entropy);
    data.volume   = bilinear_interpolator(volume);

    // Check if pK, lnK or logK was provided to calculate the standard molar Gibbs energy of the reaction
    if(not pK.empty()) data.gibbs = deltaGibssFromPK(data.pK);
    if(not lnK.empty()) data.gibbs = deltaGibssFromLnK(data.lnK);
    if(not logK.empty()) data.gibbs = deltaGibssFromLogK(data.logK);

    return data;
}

auto parseAqueousThermoDataHKF(const xml_node& node) -> boost::optional<AqueousThermoDataHKF>
{
    const bool emptyGf = node.child("Gf").text().empty();
    const bool emptyHf = node.child("Hf").text().empty();

    if(emptyGf or emptyHf)
        return boost::optional<AqueousThermoDataHKF>();

    AqueousThermoDataHKF hkf;

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

auto parseGaseousThermoDataHKF(const xml_node& node) -> boost::optional<GaseousThermoDataHKF>
{
    const bool emptyGf = node.child("Gf").text().empty();
    const bool emptyHf = node.child("Hf").text().empty();

    if(emptyGf or emptyHf)
        return boost::optional<GaseousThermoDataHKF>();

    GaseousThermoDataHKF hkf;

    hkf.Gf = node.child("Gf").text().as_double();
    hkf.Hf = node.child("Hf").text().as_double();
    hkf.Sr = node.child("Sr").text().as_double();
    hkf.a = node.child("a").text().as_double();
    hkf.b = node.child("b").text().as_double();
    hkf.c = node.child("c").text().as_double();
    hkf.Tmax = node.child("Tmax").text().as_double();

    return hkf;
}

auto parseMineralThermoDataHKF(const xml_node& node) -> boost::optional<MineralThermoDataHKF>
{
    const bool emptyGf = node.child("Gf").text().empty();
    const bool emptyHf = node.child("Hf").text().empty();

    if(emptyGf or emptyHf)
        return boost::optional<MineralThermoDataHKF>();

    MineralThermoDataHKF hkf;

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

auto parseAqueousThermoData(const xml_node& node) -> AqueousThermoData
{
    AqueousThermoData thermo_data;

    if(not node.child("interpolated").empty())
        thermo_data.interpolated = parseThermoDataInterpolated(node.child("interpolated"));

    if(not node.child("reaction").empty())
        thermo_data.reaction = parseThermoDataReaction(node.child("reaction"));

    if(not node.child("hkf").empty())
        thermo_data.hkf = parseAqueousThermoDataHKF(node.child("hkf"));

    return thermo_data;
}

auto parseGaseousThermoData(const xml_node& node) -> GaseousThermoData
{
    GaseousThermoData thermo_data;

    if(not node.child("interpolated").empty())
        thermo_data.interpolated = parseThermoDataInterpolated(node.child("interpolated"));

    if(not node.child("reaction").empty())
        thermo_data.reaction = parseThermoDataReaction(node.child("reaction"));

    if(not node.child("hkf").empty())
        thermo_data.hkf = parseGaseousThermoDataHKF(node.child("hkf"));

    return thermo_data;
}

auto parseMineralThermoData(const xml_node& node) -> MineralThermoData
{
    MineralThermoData thermo_data;

    if(not node.child("interpolated").empty())
        thermo_data.interpolated = parseThermoDataInterpolated(node.child("interpolated"));

    if(not node.child("reaction").empty())
        thermo_data.reaction = parseThermoDataReaction(node.child("reaction"));

    if(not node.child("hkf").empty())
        thermo_data.hkf = parseMineralThermoDataHKF(node.child("hkf"));

    return thermo_data;
}

auto parseSpecies(const xml_node& node, GeneralSpecies& species) -> void
{
    // Set the name of the species
    species.setName(node.child("name").text().get());

    // Set the chemical formula of the species
    species.setFormula(node.child("formula").text().get());

    // Set the elemental formula of the species
    species.setElements(node.child("elemental_formula").text().get());

    // Set the molar mass of the species
    if(not node.child("molar_mass").text().empty())
    {
        const auto value = node.child("molar_mass").text().as_double();
        const auto units = node.child("molar_mass").attribute("units").as_string();
        species.setMolarMass(units::MolarMass(value, units));
    }
}

auto parseAqueousSpecies(const xml_node& node, AqueousSpecies& species) -> void
{
    parseSpecies(node, species);

    // Set the elemental charge of the species
    species.setCharge(node.child("charge").text().as_double());

    // Parse the complex formula of the aqueous species (if any)
    species.setDissociation(node.child("dissociation").text().get());

    // Parse the thermodynamic data of the aqueous species
    species.setThermoData(parseAqueousThermoData(node.child("thermo_data")));
}

auto parseGaseousSpecies(const xml_node& node, GaseousSpecies& species) -> void
{
    parseSpecies(node, species);

    // Set the name of the gas
    species.setGas(node.child("gas").text().get());

    // Parse the thermodynamic data of the gaseous species
    species.setThermoData(parseGaseousThermoData(node.child("thermo_data")));
}

auto parseMineralSpecies(const xml_node& node, MineralSpecies& species) -> void
{
    parseSpecies(node, species);

    // Set the molar volume of the mineral species
    if(not node.child("molar_volume").text().empty())
    {
        const auto value = node.child("molar_volume").text().as_double();
        const auto units = node.child("molar_volume").attribute("units").as_string();
        species.setMolarVolume(units::MolarVolume(value, units));
    }

    // Parse the thermodynamic data of the mineral species
    species.setThermoData(parseMineralThermoData(node.child("thermo_data")));
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
        return contained(species.elementNames(), elements);
    };

    return collectSpecies(map, f);
}


} /* namespace internal */

using namespace internal;

class Database::Impl
{
private:
    /// Auxiliary types for a map of aqueous, gaseous, and mineral species
    using AqueousSpeciesMap = std::map<std::string, AqueousSpecies>;
    using GaseousSpeciesMap = std::map<std::string, GaseousSpecies>;
    using MineralSpeciesMap = std::map<std::string, MineralSpecies>;

    /// The set of all aqueous species in the database
    AqueousSpeciesMap m_aqueous_species_map;

    /// The set of all gaseous species in the database
    GaseousSpeciesMap m_gaseous_species_map;

    /// The set of all mineral species in the database
    MineralSpeciesMap m_mineral_species_map;

public:
    Impl()
    {}

    Impl(const std::string& filename)
    {
        parse(filename);
    }

    template<typename Key, typename Value>
    std::vector<Value> collectValues(const std::map<Key, Value>& map)
    {
        std::vector<Value> species;
        species.reserve(map.size());
        for(const auto& pair : map)
            species.push_back(pair.second);
        return species;
    }

    auto aqueousSpecies() -> std::vector<AqueousSpecies>
    {
        return collectValues(m_aqueous_species_map);
    }

    auto aqueousSpecies(const std::string& name) const -> const AqueousSpecies&
    {
        if(m_aqueous_species_map.count(name) == 0)
            nonexistentSpeciesError("aqueous", name);

        return m_aqueous_species_map.find(name)->second;
    }

    auto gaseousSpecies() -> std::vector<GaseousSpecies>
    {
        return collectValues(m_gaseous_species_map);
    }

    auto gaseousSpecies(const std::string& name) const -> const GaseousSpecies&
    {
        if(m_gaseous_species_map.count(name) == 0)
            nonexistentSpeciesError("gaseous", name);

        return m_gaseous_species_map.find(name)->second;
    }

    auto mineralSpecies() -> std::vector<MineralSpecies>
    {
        return collectValues(m_mineral_species_map);
    }

    auto mineralSpecies(const std::string& name) const -> const MineralSpecies&
    {
        if(m_mineral_species_map.count(name) == 0)
            nonexistentSpeciesError("mineral", name);

        return m_mineral_species_map.find(name)->second;
    }

    auto containsAqueousSpecies(const std::string& species) const -> bool
    {
        return m_aqueous_species_map.count(species) != 0;
    }

    auto containsGaseousSpecies(const std::string& species) const -> bool
    {
        return m_gaseous_species_map.count(species) != 0;
    }

    auto containsMineralSpecies(const std::string& species) const -> bool
    {
        return m_mineral_species_map.count(species) != 0;
    }

    auto aqueousSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<std::string>
    {
        return speciesWithElements(elements, m_aqueous_species_map);
    }

    auto gaseousSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<std::string>
    {
        return speciesWithElements(elements, m_gaseous_species_map);
    }

    auto mineralSpeciesWithElements(const std::vector<std::string>& elements) const -> std::vector<std::string>
    {
        return speciesWithElements(elements, m_mineral_species_map);
    }

    auto parse(const std::string& filename) -> void
    {
        xml_document doc;

        auto result = doc.load_file(filename.c_str());

        xml_node database = doc.child("database");

        if(not result)
            throw std::runtime_error("Error: could not open the database file.");

        for(xml_node node : database)
        {
            auto type = std::string(node.child("type").text().get());

            if(type == std::string("Aqueous"))
            {
                AqueousSpecies species;
                parseAqueousSpecies(node, species);

                m_aqueous_species_map[species.name()] = species;
            }
            else if(type == std::string("Gaseous"))
            {
                GaseousSpecies species;
                parseGaseousSpecies(node, species);

                m_gaseous_species_map[species.name()] = species;
            }

            else if(type == std::string("Mineral"))
            {
                MineralSpecies species;
                parseMineralSpecies(node, species);

                m_mineral_species_map[species.name()] = species;
            }
            else throw std::runtime_error("Error: the type of the species is unknown.");
        }
    }

};

Database::Database()
: pimpl(new Impl())
{}

Database::Database(const std::string& filename)
: pimpl(new Impl(filename))
{}

auto Database::aqueousSpecies() -> std::vector<AqueousSpecies>
{
    return pimpl->aqueousSpecies();
}

auto Database::aqueousSpecies(const std::string& name) const -> const AqueousSpecies&
{
    return pimpl->aqueousSpecies(name);
}

auto Database::gaseousSpecies() -> std::vector<GaseousSpecies>
{
    return pimpl->gaseousSpecies();
}

auto Database::gaseousSpecies(const std::string& name) const -> const GaseousSpecies&
{
    return pimpl->gaseousSpecies(name);
}

auto Database::mineralSpecies() -> std::vector<MineralSpecies>
{
    return pimpl->mineralSpecies();
}

auto Database::mineralSpecies(const std::string& name) const -> const MineralSpecies&
{
    return pimpl->mineralSpecies(name);
}

auto Database::containsAqueousSpecies(const std::string& species) const -> bool
{
    return pimpl->containsAqueousSpecies(species);
}

auto Database::containsGaseousSpecies(const std::string& species) const -> bool
{
    return pimpl->containsGaseousSpecies(species);
}

auto Database::containsMineralSpecies(const std::string& species) const -> bool
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

} /* namespace Reaktor */

