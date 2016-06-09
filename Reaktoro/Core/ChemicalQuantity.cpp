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

#include "ChemicalQuantity.hpp"

// C++ includes
#include <map>

// Reaktoro includes
#include <Reaktoro/Common/ConvertUtils.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalPropertiesAqueousPhase.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {
namespace {

/// The stored result of ln(10).
const double ln10 = 2.30258509299;

/// Return the default units of a quantity.
auto defaultQuantityUnits(std::string quantity) -> std::string
{
    static const std::map<std::string, std::string> default_units =
    {
        {"elementamount"        , "mol"},
        {"elementamountinphase" , "mol"},
        {"elementmass"          , "kg"},
        {"elementmassinphase"   , "kg"},
        {"elementmolality"      , "molal"},
        {"elementmolarity"      , "molar"},
        {"fluidvolume"          , "m3"},
        {"fugacity"             , "bar"},
        {"ionicstrength"        , "molal"},
        {"phaseamount"          , "mol"},
        {"phasemass"            , "kg"},
        {"phasevolume"          , "m3"},
        {"pressure"             , "pascal"},
        {"reactionrate"         , "mol/s"},
        {"solidvolume"          , "m3"},
        {"speciesamount"        , "mol"},
        {"speciesmass"          , "kg"},
        {"speciesmolality"      , "molal"},
        {"speciesmolarity"      , "molar"},
        {"temperature"          , "kelvin"},
        {"time"                 , "s"},
        {"t"                    , "s"},
        {"volume"               , "m3"},
    };
    auto iter = default_units.find(quantity);
    return iter != default_units.end() ? iter->second : "";
}

/// A type used to describe the scale of the quantity.
enum class QuantityScale
{
    one, ln, log, exp
};

/// A type used to describe the attributes collected from a formatted quantity string.
struct QuantityData
{
    /// The formated string used to create this QuantityData instance
    std::string str;

    /// The name of the quantity (e.g., temperature, elementMass, etc.)
    std::string quantity;

    /// The non-keyword arguments of the quantity function
    std::vector<std::string> args;

    /// The scale of the quantity.
    QuantityScale scale = QuantityScale::one;

    /// The units of the quantity.
    std::string units;

    /// The units conversion factor
    double factor = 1.0;
};

/// Convert a string into a QuantityScale
auto convertStringToQuantityScale(std::string scale) -> QuantityScale
{
    if(scale == "one") return QuantityScale::one;
    if(scale == "ln")  return QuantityScale::ln;
    if(scale == "log") return QuantityScale::log;
    if(scale == "exp") return QuantityScale::exp;
    RuntimeError("Could not convert the string `" +
        scale + "` into a QuantityScale.", "This is "
            "not a supported scale.")
}

/// Convert a formatted string into a QuantityData instance.
auto convertStringToQuantityData(std::string str) -> QuantityData
{
    // Remove leading and trailing white spaces
    str = trim(str);

    // Find the indices of first `(` and last `)` in the string
    auto ibracket_begin = str.find("(");
    auto ibracket_end = str.find_last_of(")");

    // Assert ( comes before )
    Assert(ibracket_begin <= ibracket_end,
        "Could not parse the given quantity string `" + str + "`.",
        "Ensure opening bracket `(` comes before closing bracket `)`.");

    // Assert ) is the final char
    if(ibracket_end < str.size())
        Assert(str.back() == ')',
            "Could not parse the given quantity string `" + str + "`.",
            "Ensure closing bracket `)` is present.");

    QuantityData data;

    // Set the formatted string
    data.str = str;

    // Extract the quantity name
    data.quantity = lowercase(str.substr(0, ibracket_begin));

    // Set the default quantity units
    data.units = defaultQuantityUnits(data.quantity);

    // Skip the rest if there are no brackets
    if(ibracket_begin == std::string::npos)
        return data;

    // Create a string with the words inside brackets
    std::string inside = str.substr(ibracket_begin + 1);
    inside.pop_back(); // remove the trailing bracket )

    // Split the inner words at space
    auto words = split(inside);

    // Loop over all inner words
    for(auto word : words)
    {
        auto pair = split(word, "=");

        if(pair.size() == 1)
        {
            data.args.push_back(word);
        }
        else if(pair.size() == 2)
        {
            if(pair[0] == "units")
            {
                if(data.quantity != "temperature")
                    data.factor = units::convert(1.0, data.units, pair[1]);
                data.units = pair[1];
            }
            else if(pair[0] == "scale")
                data.scale = convertStringToQuantityScale(pair[1]);
            else RuntimeError("Could not parse the given quantity string `" + str + "`.",
                "The provided keyword `" + pair[0] + "` is not supported.");
        }
    }

    return data;
}

/// Return the value with given scale.
auto applyQuantityScale(double val, const QuantityScale& scale) -> double
{
    switch(scale)
    {
        case QuantityScale::ln:  return std::log(val);
        case QuantityScale::log: return std::log10(val);
        case QuantityScale::exp: return std::exp(val);
        default: return val;
    }
}

} // namespace

struct ChemicalQuantity::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The reactions in the chemical system
    ReactionSystem reactions;

    /// The chemical state of the system
    ChemicalState state;

    /// The thermodynamic properties of the chemical system at (*T*, *P*, **n**)
    ChemicalProperties properties;

    /// The progress variable at which the chemical state is referred (if time, in units of s)
    double t;

    /// The temperature of the chemical system (in units of K).
    double T;

    /// The pressure of the chemical system (in units of Pa).
    double P;

    /// The molar amounts of the species in the chemical system (in units of mol).
    Vector n;

    /// The rates of the reactions in the chemical system (in units of mol/s).
    ChemicalVector r;

    /// The index of aqueous phase
    Index iAqueous;

    /// The index of aqueous species H2O
    Index iH2O;

    /// The index of aqueous species H+
    Index iH;

    /// All created chemical quantity functions from formatted strings
    std::map<std::string, Function> function_map;

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a custom Impl instance with given ChemicalSystem object
    Impl(const ChemicalSystem& system)
    : system(system)
    {
        initialize();
    }

    /// Construct a custom Impl instance with given ReactionSystem object
    Impl(const ReactionSystem& reactions)
    : system(reactions.system()), reactions(reactions)
    {
        initialize();
    }

    /// Initialize the Impl instance
    auto initialize() -> void
    {
        iAqueous = system.indexPhase("Aqueous");
        iH2O = system.indexSpeciesAny(alternativeWaterNames());
        iH = system.indexSpeciesAny(alternativeChargedSpeciesNames("H+"));
    }

    /// Update the state of the chemical quantity instance
    auto update(const ChemicalState& state) -> void
    {
        update(state, 0.0);
    }

    /// Update the state of the chemical quantity instance
    auto update(const ChemicalState& state_, double t_) -> void
    {
        // Update the chemical state of the system
        state = state_;

        // Update the progress variable
        t = t_;

        // Update the temperature, pressure and molar composition of the system
        T = state.temperature();
        P = state.pressure();
        n = state.speciesAmounts();

        // Update the thermodynamic properties of the system
        properties = system.properties(T, P, n);

        // Update the rates of the reactions
        if(!reactions.reactions().empty())
            r = reactions.rates(properties);
    }

    auto function(std::string str) -> Function
    {
        auto it = function_map.find(str);
        if(it != function_map.end())
            return it->second;

        QuantityData data = convertStringToQuantityData(str);

        Function newfunc = createFunction(data);

        function_map.insert({str, newfunc});

        return newfunc;
    }

    auto createFunction(const QuantityData& data) const -> Function
    {
        auto convert_and_apply_scale = [](double val, const QuantityData& data)
        {
            return applyQuantityScale(val * data.factor, data.scale);
        };

        auto create_function_time = [=]() -> Function
        {
            auto func = [=]() -> double
            {
                return convert_and_apply_scale(t, data);
            };
            return func;
        };

        auto create_function_temperature = [=]() -> Function
        {
            auto func = [=]() -> double
            {
                const double Tval = units::convert(T, "kelvin", data.units);
                return applyQuantityScale(Tval, data.scale);
            };
            return func;
        };

        auto create_function_pressure = [=]() -> Function
        {
            auto func = [=]() -> double
            {
                return convert_and_apply_scale(P, data);
            };
            return func;
        };

        auto create_function_elementAmount = [=]() -> Function
        {
            Assert(data.args.size() == 1,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting only one unnamed argument: the element name.");

            const Index ielement = system.indexElement(data.args[0]);

            Assert(ielement < system.numElements(),
                "Could not create the function for the quantity `" + data.str + "`.",
                "The element `" + data.args[0] + "` is not present in the system.");

            auto func = [=]() -> double
            {
                const double amount = state.elementAmount(ielement);
                return convert_and_apply_scale(amount, data);
            };
            return func;
        };

        auto create_function_elementAmountInPhase = [=]() -> Function
        {
            Assert(data.args.size() == 2,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting two unnamed arguments: the element and phase names.");

            const Index ielement = system.indexElement(data.args[0]);
            const Index iphase = system.indexPhase(data.args[1]);

            Assert(ielement < system.numElements(),
                "Could not create the function for the quantity `" + data.str + "`.",
                "The element `" + data.args[0] + "` is not present in the system.");

            Assert(iphase < system.numPhases(),
                "Could not create the function for the quantity `" + data.str + "`.",
                "The phase `" + data.args[0] + "` is not present in the system.");

            auto func = [=]() -> double
            {
                const double amount = state.elementAmountInPhase(ielement, iphase);
                return convert_and_apply_scale(amount, data);
            };
            return func;
        };

        auto create_function_elementMass = [=]() -> Function
        {
            Assert(data.args.size() == 1,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting only one unnamed argument: the element name.");

            const Index ielement = system.indexElement(data.args[0]);

            Assert(ielement < system.numElements(),
                "Could not create the function for the quantity `" + data.str + "`.",
                "The element `" + data.args[0] + "` is not present in the system.");

            const double molar_mass = system.element(ielement).molarMass();

            auto func = [=]() -> double
            {
                const double amount = state.elementAmount(ielement);
                const double mass = amount * molar_mass;
                return convert_and_apply_scale(mass, data);
            };
            return func;
        };

        auto create_function_elementMassInPhase = [=]() -> Function
        {
            Assert(data.args.size() == 2,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting two unnamed arguments: the element and phase names.");

            const Index ielement = system.indexElement(data.args[0]);
            const Index iphase = system.indexPhase(data.args[1]);

            Assert(ielement < system.numElements(),
                "Could not create the function for the quantity `" + data.str + "`.",
                "The element `" + data.args[0] + "` is not present in the system.");

            Assert(iphase < system.numPhases(),
                "Could not create the function for the quantity `" + data.str + "`.",
                "The phase `" + data.args[0] + "` is not present in the system.");

            const double molar_mass = system.element(ielement).molarMass();

            auto func = [=]() -> double
            {
                const double amount = state.elementAmountInPhase(ielement, iphase);
                const double mass = amount * molar_mass;
                return convert_and_apply_scale(mass, data);
            };
            return func;
        };

        auto create_function_elementMolality = [=]() -> Function
        {
            Assert(data.args.size() == 1,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting only one unnamed argument: the element name.");

            const Index ielement = system.indexElement(data.args[0]);

            Assert(ielement < system.numElements(),
                "Could not create the function for the quantity `" + data.str + "`.",
                "The element `" + data.args[0] + "` is not present in the system.");

            // Return zero if no aqueous phase or water species
            if(iAqueous >= system.numPhases() || iH2O >= system.numSpecies())
                return []() { return 0.0; };

            auto func = [=]() -> double
            {
                const double amount = state.elementAmountInPhase(ielement, iAqueous);
                const double kgH2O = state.speciesAmount(iH2O) * waterMolarMass;
                const double mi = kgH2O ? amount/kgH2O : 0.0;
                return convert_and_apply_scale(mi, data);
            };
            return func;
        };

        auto create_function_elementMolarity = [=]() -> Function
        {
            Assert(data.args.size() == 1,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting only one unnamed argument: the element name.");

            const Index ielement = system.indexElement(data.args[0]);

            Assert(ielement < system.numElements(),
                "Could not create the function for the quantity `" + data.str + "`.",
                "The element `" + data.args[0] + "` is not present in the system.");

            // Return zero if no aqueous phase or water species
            if(iAqueous >= system.numPhases() || iH2O >= system.numSpecies())
                return []() { return 0.0; };

            auto func = [=]() -> double
            {
                const double amount = state.elementAmountInPhase(ielement, iAqueous);
                const double volume = properties.phaseVolumes()[iAqueous].val;
                const double liter = convertCubicMeterToLiter(volume);
                const double ci = liter ? amount/liter : 0.0;
                return convert_and_apply_scale(ci, data);
            };
            return func;
        };

        auto create_function_speciesAmount = [=]() -> Function
        {
            Assert(data.args.size() == 1,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting only one unnamed argument: the species name.");

            const Index ispecies = system.indexSpecies(data.args[0]);

            Assert(ispecies < system.numSpecies(),
                "Could not create the function for the quantity `" + data.str + "`.",
                "The species `" + data.args[0] + "` is not present in the system.");

            auto func = [=]() -> double
            {
                const double amount = state.speciesAmount(ispecies);
                return convert_and_apply_scale(amount, data);
            };
            return func;
        };

        auto create_function_speciesMass = [=]() -> Function
        {
            Assert(data.args.size() == 1,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting only one unnamed argument: the species name.");

            const Index ispecies = system.indexSpecies(data.args[0]);

            Assert(ispecies < system.numSpecies(),
                "Could not create the function for the quantity `" + data.str + "`.",
                "The species `" + data.args[0] + "` is not present in the system.");

            const double molar_mass = system.species(ispecies).molarMass();

            auto func = [=]() -> double
            {
                const double amount = state.speciesAmount(ispecies);
                const double mass = amount * molar_mass;
                return convert_and_apply_scale(mass, data);
            };
            return func;
        };

        auto create_function_speciesMolality = [=]() -> Function
        {
            Assert(data.args.size() == 1,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting only one unnamed argument: the species name.");

            const Index ispecies = system.indexSpecies(data.args[0]);

            Assert(ispecies < system.numSpecies(),
                "Could not create the function for the quantity `" + data.str + "`.",
                "The species `" + data.args[0] + "` is not present in the system.");

            // Return zero if no aqueous phase or water species
            if(iAqueous >= system.numPhases() || iH2O >= system.numSpecies())
                return []() { return 0.0; };

            const Index index_first_aqueous_species = system.indexFirstSpeciesInPhase(iAqueous);
            const Index num_aqueous_species = system.numSpeciesInPhase(iAqueous);

            Assert(ispecies - index_first_aqueous_species < num_aqueous_species,
                "Could not create the function for the quantity `" + data.str + "`.",
                "The species `" + data.args[0] + "` is not present in the aqueous phase.");

            auto func = [=]() -> double
            {
                const double amount = state.speciesAmount(ispecies);
                const double kgH2O = state.speciesAmount(iH2O) * waterMolarMass;
                const double mi = kgH2O ? amount/kgH2O : 0.0;
                return convert_and_apply_scale(mi, data);
            };
            return func;
        };

        auto create_function_speciesMolarity = [=]() -> Function
        {
            Assert(data.args.size() == 1,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting only one unnamed argument: the species name.");

            const Index ispecies = system.indexSpecies(data.args[0]);

            Assert(ispecies < system.numSpecies(),
                "Could not create the function for the quantity `" + data.str + "`.",
                "The species `" + data.args[0] + "` is not present in the system.");

            // Return zero if no aqueous phase or water species
            if(iAqueous >= system.numPhases() || iH2O >= system.numSpecies())
                return []() { return 0.0; };

            const Index index_first_aqueous_species = system.indexFirstSpeciesInPhase(iAqueous);
            const Index num_aqueous_species = system.numSpeciesInPhase(iAqueous);

            Assert(ispecies - index_first_aqueous_species < num_aqueous_species,
                "Could not create the function for the quantity `" + data.str + "`.",
                "The species `" + data.args[0] + "` is not present in the aqueous phase.");

            auto func = [=]() -> double
            {
                const double amount = state.speciesAmount(ispecies);
                const double volume = properties.phaseVolumes()[iAqueous].val;
                const double liter = convertCubicMeterToLiter(volume);
                const double ci = liter ? amount/liter : 0.0;
                return convert_and_apply_scale(ci, data);
            };
            return func;
        };

        auto create_function_speciesMolarFraction = [=]() -> Function
        {
            Assert(data.args.size() == 1,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting only one unnamed argument: the species name.");

            const Index ispecies = system.indexSpecies(data.args[0]);

            Assert(ispecies < system.numSpecies(),
                "Could not create the function for the quantity `" + data.str + "`.",
                "The species `" + data.args[0] + "` is not present in the system.");

            auto func = [=]() -> double
            {
                const double xi = properties.molarFractions().val[ispecies];
                return convert_and_apply_scale(xi, data);
            };
            return func;
        };

        auto create_function_speciesActivity = [=]() -> Function
        {
            Assert(data.args.size() == 1,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting only one unnamed argument: the species name.");

            const Index ispecies = system.indexSpecies(data.args[0]);

            Assert(ispecies < system.numSpecies(),
                "Could not create the function for the quantity `" + data.str + "`.",
                "The species `" + data.args[0] + "` is not present in the system.");

            auto func = [=]() -> double
            {
                const double ln_ai = properties.lnActivities().val[ispecies];
                return convert_and_apply_scale(std::exp(ln_ai), data);
            };
            return func;
        };

        auto create_function_speciesActivityCoefficient = [=]() -> Function
        {
            Assert(data.args.size() == 1,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting only one unnamed argument: the species name.");

            const Index ispecies = system.indexSpecies(data.args[0]);

            Assert(ispecies < system.numSpecies(),
                "Could not create the function for the quantity `" + data.str + "`.",
                "The species `" + data.args[0] + "` is not present in the system.");

            auto func = [=]() -> double
            {
                const double ln_gi = properties.lnActivityCoefficients().val[ispecies];
                return convert_and_apply_scale(std::exp(ln_gi), data);
            };
            return func;
        };

        auto create_function_phaseAmount = [=]() -> Function
        {
            Assert(data.args.size() == 1,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting only one unnamed argument: the phase name.");

            const Index iphase = system.indexPhase(data.args[0]);

            Assert(iphase < system.numPhases(),
                "Could not create the function for the quantity `" + data.str + "`.",
                "The phase `" + data.args[0] + "` is not present in the system.");

            auto func = [=]() -> double
            {
                const double amount = state.phaseAmount(iphase);
                return convert_and_apply_scale(amount, data);
            };
            return func;
        };

        auto create_function_phaseMass = [=]() -> Function
        {
            Assert(data.args.size() == 1,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting only one unnamed argument: the phase name.");

            const Index iphase = system.indexPhase(data.args[0]);

            Assert(iphase < system.numPhases(),
                "Could not create the function for the quantity `" + data.str + "`.",
                "The phase `" + data.args[0] + "` is not present in the system.");

            auto func = [=]() -> double
            {
                const double mass = properties.phaseMasses().val[iphase];
                return convert_and_apply_scale(mass, data);
            };
            return func;
        };

        auto create_function_phaseVolume = [=]() -> Function
        {
            Assert(data.args.size() == 1,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting only one unnamed argument: the phase name.");

            const Index iphase = system.indexPhase(data.args[0]);

            Assert(iphase < system.numPhases(),
                "Could not create the function for the quantity `" + data.str + "`.",
                "The phase `" + data.args[0] + "` is not present in the system.");

            auto func = [=]() -> double
            {
                const double volume = properties.phaseVolumes().val[iphase];
                return convert_and_apply_scale(volume, data);
            };
            return func;
        };

        auto create_function_fluidVolume = [=]() -> Function
        {
            Assert(data.args.size() == 0,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting no unnamed arguments.");

            auto func = [=]() -> double
            {
                const double volume = properties.fluidVolume().val;
                return convert_and_apply_scale(volume, data);
            };
            return func;
        };

        auto create_function_solidVolume = [=]() -> Function
        {
            Assert(data.args.size() == 0,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting no unnamed arguments.");

            auto func = [=]() -> double
            {
                const double volume = properties.solidVolume().val;
                return convert_and_apply_scale(volume, data);
            };
            return func;
        };

        auto create_function_porosity = [=]() -> Function
        {
            Assert(data.args.size() == 0,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting no unnamed arguments.");

            auto func = [=]() -> double
            {
                const double totalVolume = properties.volume().val;
                const double solidVolume = properties.solidVolume().val;
                const double porosity = 1 - solidVolume/totalVolume;
                return applyQuantityScale(porosity, data.scale);
            };
            return func;
        };

        auto create_function_ionicStrength = [=]() -> Function
        {
            Assert(data.args.size() == 0,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting no unnamed arguments.");

            // Return zero if there are no aqueous phase
            if(iAqueous >= system.numPhases())
                return []() { return 0.0; };

            auto func = [=]() -> double
            {
                const double I = properties.aqueous().ionicStrength().val;
                return convert_and_apply_scale(I, data);
            };
            return func;
        };

        auto create_function_pH = [=]() -> Function
        {
            Assert(data.args.size() == 0,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting no unnamed arguments.");

            // Return zero if no aqueous phase or hydron species
            if(iAqueous >= system.numPhases() || iH >= system.numSpecies())
                return []() { return 0.0; };

            auto func = [=]() -> double
            {
                const double ln_aH = properties.lnActivities().val[iH];
                const double pH = -ln_aH/ln10;
                return applyQuantityScale(pH, data.scale);
            };
            return func;
        };

        auto create_function_pE = [=]() -> Function
        {
            Assert(data.args.size() == 0,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting no unnamed arguments.");

            // Return zero if there are no aqueous phase
            if(iAqueous >= system.numPhases())
                return []() { return 0.0; };

            auto func = [=]() -> double
            {
                const double pE = properties.aqueous().pE().val;
                return applyQuantityScale(pE, data.scale);
            };
            return func;
        };

        auto create_function_Eh = [=]() -> Function
        {
            Assert(data.args.size() == 0,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting no unnamed arguments.");

            // Return zero if there are no aqueous phase
            if(iAqueous >= system.numPhases())
                return []() { return 0.0; };

            auto func = [=]() -> double
            {
                const double Eh = properties.aqueous().Eh().val;
                return applyQuantityScale(Eh, data.scale);
            };
            return func;
        };

        auto create_function_reactionRate = [=]() -> Function
        {
            Assert(data.args.size() == 1,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting only one unnamed argument: the reaction name.");

            const Index ireaction = reactions.indexReaction(data.args[0]);

            Assert(ireaction < reactions.numReactions(),
                "Could not create the function for the quantity `" + data.str + "`.",
                "The reaction `" + data.args[0] + "` is not present in the reaction system.");

            auto func = [=]() -> double
            {
                const double ri = r.val[ireaction];
                return convert_and_apply_scale(ri, data);
            };
            return func;
        };

        auto create_function_reactionEquilibriumIndex = [=]() -> Function
        {
            Assert(data.args.size() == 1,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting only one unnamed argument: the reaction name.");

            const Index ireaction = reactions.indexReaction(data.args[0]);

            Assert(ireaction < reactions.numReactions(),
                "Could not create the function for the quantity `" + data.str + "`.",
                "The reaction `" + data.args[0] + "` is not present in the reaction system.");

            auto func = [=]() -> double
            {
                const double ln_omega = reactions.reaction(ireaction).lnEquilibriumIndex(properties).val;
                return convert_and_apply_scale(std::exp(ln_omega), data);
            };
            return func;
        };

        // Create a dictionary of quantity-function pairs
        std::map<std::string, std::function<Function()>> create_function_map =
        {
            {"t"                          , create_function_time},
            {"time"                       , create_function_time},
            {"progress"                   , create_function_time},
            {"temperature"                , create_function_temperature},
            {"pressure"                   , create_function_pressure},
            {"elementamount"              , create_function_elementAmount},
            {"elementamountinphase"       , create_function_elementAmountInPhase},
            {"elementmass"                , create_function_elementMass},
            {"elementmassinphase"         , create_function_elementMassInPhase},
            {"elementmolality"            , create_function_elementMolality},
            {"elementmolarity"            , create_function_elementMolarity},
            {"speciesamount"              , create_function_speciesAmount},
            {"speciesmass"                , create_function_speciesMass},
            {"speciesmolality"            , create_function_speciesMolality},
            {"speciesmolarity"            , create_function_speciesMolarity},
            {"speciesmolarfraction"       , create_function_speciesMolarFraction},
            {"speciesactivity"            , create_function_speciesActivity},
            {"speciesactivitycoefficient" , create_function_speciesActivityCoefficient},
            {"phaseamount"                , create_function_phaseAmount},
            {"phasemass"                  , create_function_phaseMass},
            {"phasevolume"                , create_function_phaseVolume},
            {"fluidvolume"                , create_function_fluidVolume},
            {"solidvolume"                , create_function_solidVolume},
            {"porosity"                   , create_function_porosity},
            {"ionicstrength"              , create_function_ionicStrength},
            {"ph"                         , create_function_pH},
            {"pE"                         , create_function_pE},
            {"Eh"                         , create_function_Eh},
            {"reactionrate"               , create_function_reactionRate},
            {"reactionequilibriumindex"   , create_function_reactionEquilibriumIndex},
        };

        std::string quantity = lowercase(data.quantity);

        Assert(create_function_map.count(quantity),
            "Could not create the function for the quantity `" + data.str + "`.",
            "The quantity name is not supported.");

        return create_function_map.at(quantity)();
    }

    auto value(std::string str) -> double
    {
        return function(str)();
    }
};

ChemicalQuantity::ChemicalQuantity()
: pimpl(new Impl())
{}

ChemicalQuantity::ChemicalQuantity(const ChemicalQuantity& other)
: pimpl(new Impl(*other.pimpl))
{}

ChemicalQuantity::ChemicalQuantity(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

ChemicalQuantity::ChemicalQuantity(const ReactionSystem& reactions)
: pimpl(new Impl(reactions))
{}

ChemicalQuantity::~ChemicalQuantity()
{}

auto ChemicalQuantity::operator=(ChemicalQuantity other) -> ChemicalQuantity&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto ChemicalQuantity::update(const ChemicalState& state) -> void
{
    pimpl->update(state);
}

auto ChemicalQuantity::update(const ChemicalState& state, double t) -> void
{
    pimpl->update(state, t);
}

auto ChemicalQuantity::value(std::string str) const -> double
{
    return pimpl->value(str);
}

auto ChemicalQuantity::function(std::string str) const -> Function
{
    return pimpl->function(str);
}

auto ChemicalQuantity::operator[](std::string quantity) const -> double
{
    return value(quantity);
}

} // namespace Reaktoro
