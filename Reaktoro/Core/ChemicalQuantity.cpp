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
        {"elementamount"        , " mol"},
        {"elementamountinphase" , " mol"},
        {"elementmass"          , " kg"},
        {"elementmassinphase"   , " kg"},
        {"elementmolality"      , " molal"},
        {"elementmolarity"      , " molar"},
        {"fluidvolume"          , " m3"},
        {"fugacity"             , " bar"},
        {"phaseamount"          , " mol"},
        {"phasemass"            , " kg"},
        {"pressure"             , " pascal"},
        {"reactionrate"         , " mol/s"},
        {"solidvolume"          , " m3"},
        {"speciesamount"        , " mol"},
        {"speciesmass"          , " kg"},
        {"speciesmolality"      , " molal"},
        {"speciesmolarity"      , " molar"},
        {"temperature"          , " kelvin"},
        {"time"                 , " s"},
        {"volume"               , " m3"},
    };
    auto iter = default_units.find(quantity);
    return iter != default_units.end() ? iter->second : "";
}

/// A type used to describe the chemical quantity type.
enum class QuantityType
{
    Activity,
    ActivityCoefficient,
    Amount,
    Eh,
    EquilibriumIndex,
    Fugacity,
    Mass,
    Molality,
    MolarFraction,
    Molarity,
    pe,
    pH,
    Pressure,
    Progress,
    Rate,
    Temperature,
    Time,
    NotSupported
};

/// A type used to describe operators to be applied to the chemical quantity.
enum class OperatorType
{
    None, ln, log, exp
};

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
        Assert(str.back() == ")",
            "Could not parse the given quantity string `" + str + "`.",
            "Ensure closing bracket `)` is present.");

    QuantityData data;

    // Set the formatted string
    data.str = str;

    // Extract the quantity name
    data.quantity = lowercase(str.substr(0, ibracket_begin));

    // Set the default quantity units
    data.units = defaultQuantityUnits(data.quantity);

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
auto applyQuantityScale(double val, const QuantityScale& scale) const -> double
{
    switch(scale)
    {
        case OperatorType::ln:  return std::log(val);
        case OperatorType::log: return std::log10(val);
        case OperatorType::exp: return std::exp(val);
        default: return val;
    }
}

/// A type used to describe a chemical quantity.
struct Description
{
    /// The type of the quantity.
    QuantityType quantity = QuantityType::NotSupported;

    /// The index of the element for which the quantity is related.
    /// This can be empty if the quantity is not related to an element.
    Index element;

    /// The index of the phase with the element for which the quantity is related.
    /// This can be empty if the quantity is not related to an element in a phase.
    Index phase_with_element;

    /// The index of the species for which the quantity is related.
    /// This can be empty if the quantity is not related to a species.
    Index species;

    /// The index of the phase for which the quantity is related.
    /// This can be empty if the quantity is not related to a phase.
    Index phase;

    /// The index of the reaction for which the quantity is related.
    /// This can be empty if the quantity is not related to a reaction.
    Index reaction;

    /// The factor used to convert default units to desired units.
    double factor = 1.0;

    /// The temperature units if this quantity is temperature
    std::string tunits;

    /// The type of the operator to be applied to the quantity.
    OperatorType op;
};

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

    auto parse(std::string str) const -> Description
    {
        // Split the string (at spaces) into words
        auto words = split(str);

        // Check the string is not empty or contains only white spaces
        Assert(words.size() > 1,
            "Could not parse the given quantity string `" + str + "`.",
            "The string is either empty or contains only spaces.");

        // Auxiliary strings that are components of the formatted string `str`
        Description description;

        // Set the quantity name as the first word
        std::string quantity;
        quantity = words[0];

        // For each subsequent words, we expect them as `keyword=value`
        for(auto word : words)
        {
            // Split at equal sign `=`
            auto subwords = split(word, "=");

            // Assert there are two subwords
            Assert(subwords.size() == 2,
                "Could not parse the given quantity string `" + str + "`.",
                "Ensure the equal sign `=` is used correctly as `keyword=value`.");

            if(subwords[0] == "element")
            {
                // Split on `:` to check if an element in a phase is given
                auto pair = split(subwords[1], ":");

                // Check cases such as Aqueous:Ca
                if(pair.size() == 2)
                {
                    description.element = system.indexElementWithError(pair[1]);
                    description.phase_with_element = system.indexPhaseWithError(pair[0]);
                }
                // Check if only one element name is given
                else if(pair.size() == 1)
                {
                    description.element = system.indexElementWithError(subwords[1]);
                }
                // Check unsupported cases such as Aqueous:Ca:Mg
                else RuntimeError("Could not parse the given quantity string `" + str + "`.",
                    "Ensure the correct notation `phase:element` for specifying an element "
                    "in a phase (e.g., Aqueous:Ca, Gaseous:C).");
            }
            else if(subwords[0] == "species")
                description.species = system.indexSpeciesWithError(subwords[1]);
            else if(subwords[0] == "phase")
                description.phase = system.indexPhaseWithError(subwords[1]);
            else if(subwords[0] == "reaction")
                description.reaction = reactions.indexReactionWithError(subwords[1]);
            else if(subwords[0] == "units")
                if(quantity == "temperature")
                    description.tunits = subwords[1];
                else
                    description.factor = units::convert(1.0, subwords[1], defaultQuantityUnits(quantity));
            else if(subwords[0] == "scale")
                description.op = subwords[1];
            else RuntimeError("Could not parse the given quantity string `" + str + "`.",
                "The provided keyword `" + subwords[0] + "` is not supported.");
        }

        // Extract the scale of the quantity value, if any, and remove the first word from words
        if(words[0] == "ln" || words[0] == "log" || words[0] == "exp")
        {
            scale = words[0];
            words.erase(words.begin());
        }

        // Check if the string has a quantity name
        Assert(words.size() > 1,
            "Could not parse the given quantity string `" + str + "`.",
            "The string does not contain a quantity name.");

        // Set the quantity name and remove the first word from words
        quantity = words[0];
        words.erase(words.begin());

        auto create_elementAmount = [&]()
        {

        };

        // Extract the units (or non-units quantity name) from the formatted string
        auto pos = str.find("(");
        units = str.substr(0, pos);

        // Extract the entity, if any, from the formatted string
        if(pos != std::string::npos)
            entity = str.substr(pos+1, str.size()-pos-2);

        // Check if entity is in fact of type `element'phase`, e.g., `C'Aqueous`
        if(entity.size())
        {
            auto words = split(entity, "'");
            entity = words[0];
            phase  = words.size() == 2 ? words[1] : "";
        }

        //-----------------------------------------------------------------------------------

        // The map from base units to quantity keyword
        static const std::map<std::string, QuantityType> quantities_with_units =
        {
            {"s"      , QuantityType::Time},
            {"kelvin" , QuantityType::Temperature},
            {"pascal" , QuantityType::Pressure},
            {"mol"    , QuantityType::Amount},
            {"kg"     , QuantityType::Mass},
            {"molal"  , QuantityType::Molality},
            {"molar"  , QuantityType::Molarity},
            {"mol/s"  , QuantityType::Rate},
        };

        // The map of quantities that have no units
        static const std::map<std::string, QuantityType> quantities_without_units =
        {
            {"activity"            , QuantityType::Activity},
            {"activitycoefficient" , QuantityType::ActivityCoefficient},
            {"eh"                  , QuantityType::Eh},
            {"equilibriumindex"    , QuantityType::EquilibriumIndex},
            {"fugacity"            , QuantityType::Fugacity},
            {"molarfraction"       , QuantityType::MolarFraction},
            {"pe"                  , QuantityType::pe},
            {"ph"                  , QuantityType::pH},
            {"t"                   , QuantityType::Progress},
        };

        // The map of operators
        static const std::map<std::string, OperatorType> operators =
        {
            {""    , OperatorType::None},
            {"ln"  , OperatorType::ln},
            {"log" , OperatorType::log},
            {"exp" , OperatorType::exp},
        };

        Description desc;

        // Set the mathematical operator
        if(operators.count(scale))
            desc.op = operators.at(scale);
        else RuntimeError("Could not identify operator type from `" + scale + "`.",
            "This is not a valid mathematical operator.");

        // Check if the quantity has no units
        std::string lunits = lowercase(units);
        if(quantities_without_units.count(lunits))
            desc.quantity = quantities_without_units.at(lunits);

        // Check if there is an entry in the map whose key is the given units
        else if(quantities_with_units.count(units))
            desc.quantity = quantities_with_units.at(units);

        // Return the quantity keyword corresponding to a unit that is convertible to the given one
        else for(const auto& pair : quantities_with_units)
            if(units::convertible(units, pair.first)) {
                desc.quantity = pair.second;
                if(desc.quantity == QuantityType::Temperature)
                    desc.tunits = units;
                else
                    desc.factor = units::convert(1.0, pair.first, units);
                break;
            }

        // Error identifying the quantity
        Assert(desc.quantity != QuantityType::NotSupported,
            "Could not identify the quantity type from `" + units + "`.",
            "This is neither a valid quantity name nor a supported quantity units.");

        // Set the appropriate entity name
        desc.element  = system.indexElement(entity);
        desc.species  = system.indexSpecies(entity);
        desc.phase    = system.indexPhase(entity);
        desc.reaction = reactions.indexReaction(entity);

        // Overwrite the index of phase if a phase is provided
        if(!phase.empty())
            desc.phase = system.indexPhase(phase);

        //-----------------------------------------------------------------------------------

        // Validate the formatted string
        switch(desc.quantity)
        {
        case QuantityType::pH:
        case QuantityType::pe:
        case QuantityType::Eh:
        case QuantityType::Time:
        case QuantityType::Pressure:
        case QuantityType::Temperature:
        case QuantityType::Progress:
            break;

        case QuantityType::Mass:
        case QuantityType::Amount:
            Assert(desc.species < system.numSpecies() ||
                   desc.element < system.numElements() ||
                    desc.phase  < system.numPhases(),
                    "Could not parse the quantity string `" + str + "` with entity `" + entity + "`.",
                    "The quantity `" + units + "` requires an existing "
                        "element, species, or phase name, e.g., `mol(H)`, `kg(Calcite)`.");
            break;

        case QuantityType::Molality:
        case QuantityType::Molarity:
            Assert(desc.species < system.numSpecies() ||
                   desc.element < system.numElements(),
                    "Could not parse the quantity string `" + str + "` with entity `" + entity + "`.",
                    "The quantity `" + units + "` requires an existing "
                        "element or species name, e.g., `molal(C)`, `molar(Ca++)`.");
            break;

        case QuantityType::MolarFraction:
        case QuantityType::Activity:
        case QuantityType::ActivityCoefficient:
        case QuantityType::Fugacity:
            Assert(desc.species < system.numSpecies(),
                "Could not parse the quantity string `" + str + "` "
                    "with species name `" + entity + "`.",
                        "There is no species with name `" + entity + "`.");
            break;

        case QuantityType::Rate:
        case QuantityType::EquilibriumIndex:
            Assert(desc.reaction < reactions.numReactions(),
                "Could not parse the quantity string `" + str + "` "
                    "with reaction name `" + entity + "`.",
                        "There is no reaction with name `" + entity + "`.");
            break;
        default:
            break;
        }

        return desc;
    }

    auto apply(double val, const Description& desc) const -> double
    {
        if(desc.quantity == QuantityType::Temperature)
            val = units::convert(val, "kelvin", desc.tunits);
        else val *= desc.factor;

        switch(desc.op) {
            case OperatorType::ln:  return std::log(val);
            case OperatorType::log: return std::log10(val);
            case OperatorType::exp: return std::exp(val);
            default: return val; }
    }

    auto function(std::string str) -> Function
    {
        auto it = function_map.find(str);
        if(it != function_map.end())
            return it->second;

        Function newfunc = functionaux(parse(str));

        function_map.insert({str, newfunc});

        return newfunc;
    }

    auto functionaux(const Description& desc) const -> Function
    {
        auto time = [=]() -> double
        {
            return apply(t, desc);
        };

        auto progress = [=]() -> double
        {
            return apply(t, desc);
        };

        auto temperature = [=]() -> double
        {
            return apply(T, desc);
        };

        auto pressure = [=]() -> double
        {
            return apply(P, desc);
        };

        auto element_amount = [=]() -> double
        {
            return apply(state.elementAmount(desc.element), desc);
        };

        auto element_mass = [=]() -> double
        {
            const double molar_mass = system.element(desc.element).molarMass();
            const double amount = state.elementAmount(desc.element);
            return apply(amount * molar_mass, desc);
        };

        auto element_amount_in_phase = [=]() -> double
        {
            return apply(state.elementAmountInPhase(desc.element, desc.phase), desc);
        };

        auto element_mass_in_phase = [=]() -> double
        {
            const double molar_mass = system.element(desc.element).molarMass();
            const double amount = state.elementAmountInPhase(desc.element, desc.phase);
            return apply(amount * molar_mass, desc);
        };

        auto element_molality = [=]() -> double
        {
            // Return zero if no water species
            if(iH2O >= system.numSpecies()) return 0.0;
            const double amount = state.elementAmountInPhase(desc.element, iAqueous);
            const double kgH2O = state.speciesAmount(iH2O) * waterMolarMass;
            const double mi = kgH2O ? amount/kgH2O : 0.0;
            return apply(mi, desc);
        };

        auto element_molarity = [=]() -> double
        {
            // Return zero if no aqueous phase
            if(iAqueous >= system.numPhases()) return 0.0;
            const double amount = state.elementAmountInPhase(desc.element, iAqueous);
            const double volume = properties.phaseVolumes()[iAqueous].val;
            const double liter = convertCubicMeterToLiter(volume);
            const double ci = liter ? amount/liter : 0.0;
            return apply(ci, desc);
        };

        auto species_amount = [=]() -> double
        {
            return apply(state.speciesAmount(desc.species), desc);
        };

        auto species_mass = [=]() -> double
        {
            const double molar_mass = system.species(desc.species).molarMass();
            const double amount = state.speciesAmount(desc.species);
            return apply(amount * molar_mass, desc);
        };

        auto species_molality = [=]() -> double
        {
            // Return zero if no water species
            if(iH2O >= system.numSpecies()) return 0.0;
            const double amount = state.speciesAmount(desc.species);
            const double kgH2O = state.speciesAmount(iH2O) * waterMolarMass;
            const double mi = kgH2O ? amount/kgH2O : 0.0;
            return apply(mi, desc);
        };

        auto species_molarity = [=]() -> double
        {
            // Return zero if no aqueous phase
            if(iAqueous >= system.numPhases()) return 0.0;
            const double amount = state.speciesAmount(desc.species);
            const double volume = properties.phaseVolumes()[iAqueous].val;
            const double liter = convertCubicMeterToLiter(volume);
            const double ci = liter ? amount/liter : 0.0;
            return apply(ci, desc);
        };

        auto species_molarfraction = [=]() -> double
        {
            const double xi = properties.molarFractions().val[desc.species];
            return apply(xi, desc);
        };

        auto species_activity = [=]() -> double
        {
            const double ln_ai = properties.lnActivities().val[desc.species];
            return apply(std::exp(ln_ai), desc);
        };

        auto species_activitycoefficient = [=]() -> double
        {
            const double ln_gi = properties.lnActivityCoefficients().val[desc.species];
            return apply(std::exp(ln_gi), desc);
        };

        auto phase_amount = [=]() -> double
        {
            return apply(state.phaseAmount(desc.phase), desc);
        };

        auto phase_mass = [=]() -> double
        {
            const double mass = properties.phaseMasses().val[desc.phase];
            return apply(mass, desc);
        };

        auto ph = [=]() -> double
        {
            // Return zero if no hydron species
            if(iH >= system.numSpecies()) return 0.0;
            const double ln_aH = properties.lnActivities().val[iH];
            return -ln_aH/ln10;
        };

        auto reaction_rate = [=]() -> double
        {
            // Return zero if there are no reactions
            if(r.val.rows() == 0) return 0.0;
            const double ri = r.val[desc.reaction];
            return apply(ri, desc);
        };

        auto reaction_equilibriumindex = [=]() -> double
        {
            // Return zero if there are no reactions
            if(r.val.rows() == 0) return 0.0;
            const double ln_omega = reactions.reaction(desc.reaction).lnEquilibriumIndex(properties).val;
            return apply(std::exp(ln_omega), desc);
        };

        switch(desc.quantity)
        {
        case QuantityType::Time: return time;
        case QuantityType::Progress: return progress;
        case QuantityType::Temperature: return temperature;
        case QuantityType::Pressure: return pressure;
        case QuantityType::Amount:
            if(desc.element < system.numElements() && desc.phase < system.numPhases())
                return element_amount_in_phase;
            if(desc.element < system.numElements())
                return element_amount;
            if(desc.species < system.numSpecies())
                return species_amount;
            if(desc.phase < system.numPhases())
                return phase_amount;
            break;
        case QuantityType::Mass:
            if(desc.element < system.numElements() && desc.phase < system.numPhases())
                return element_mass_in_phase;
            if(desc.element < system.numElements())
                return element_mass;
            if(desc.species < system.numSpecies())
                return species_mass;
            if(desc.phase < system.numPhases())
                return phase_mass;
            break;
        case QuantityType::Molality:
            if(desc.element < system.numElements())
                return element_molality;
            if(desc.species < system.numSpecies())
                return species_molality;
            break;
        case QuantityType::Molarity:
            if(desc.element < system.numElements())
                return element_molarity;
            if(desc.species < system.numSpecies())
                return species_molarity;
            break;
        case QuantityType::Activity: return species_activity;
        case QuantityType::ActivityCoefficient: return species_activitycoefficient;
        case QuantityType::MolarFraction: return species_molarfraction;
        case QuantityType::Fugacity: return species_activity;
        case QuantityType::pH: return ph;
        case QuantityType::Rate: return reaction_rate;
        case QuantityType::EquilibriumIndex: return reaction_equilibriumindex;
        default: break;
        }
        return []() { return 0.0; };
    }

    auto createFunction(const QuantityData& data) const -> Function
    {
        auto convert_and_apply_scale = [](double val, const QuantityData& data)
        {
            return applyQuantityScale(val * data.factor, data.scale);
        };

        auto create_function_time = [=]()
        {
            auto func = [=]() -> double
            {
                return convert_and_apply_scale(t, data);
            };
            return func;
        };

        auto create_function_temperature = [=]()
        {
            auto func = [=]() -> double
            {
                const double Tval = units::convert(T, "kelvin", data.units);
                return applyQuantityScale(Tval, data.scale);
            };
            return func;
        };

        auto create_function_pressure = [=]()
        {
            auto func = [=]() -> double
            {
                return convert_and_apply_scale(P, data);
            };
            return func;
        };

        auto create_function_elementAmount = [=]() -> double
        {
            Assert(data.args.size() == 1,
                "Could not create the function for the quantity `" + data.str + "`.",
                "Expecting only one unnamed argument with the element name.");

            const Index ielement = system.indexElement(data.args[0]);

            Assert(ielement < system.numElements(),
                "Could not create the function for the quantity `" + data.str + "`.",
                "The element name `" + data.args[0] + "` is not present in the system.");

            auto func = [=]() -> double
            {
                const double bval = state.elementAmount(data.args[0]);
                return convert_and_apply_scale(bval, data);
            };
            return func;

        };

        auto element_mass = [=]() -> double
        {
            const double molar_mass = system.element(data.element).molarMass();
            const double amount = state.elementAmount(data.element);
            return convert_and_apply_scale(amount * molar_mass, data);
        };

        auto element_amount_in_phase = [=]() -> double
        {
            return convert_and_apply_scale(state.elementAmountInPhase(data.element, data.phase), data);
        };

        auto element_mass_in_phase = [=]() -> double
        {
            const double molar_mass = system.element(data.element).molarMass();
            const double amount = state.elementAmountInPhase(data.element, data.phase);
            return convert_and_apply_scale(amount * molar_mass, data);
        };

        auto element_molality = [=]() -> double
        {
            // Return zero if no water species
            if(iH2O >= system.numSpecies()) return 0.0;
            const double amount = state.elementAmountInPhase(data.element, iAqueous);
            const double kgH2O = state.speciesAmount(iH2O) * waterMolarMass;
            const double mi = kgH2O ? amount/kgH2O : 0.0;
            return convert_and_apply_scale(mi, data);
        };

        auto element_molarity = [=]() -> double
        {
            // Return zero if no aqueous phase
            if(iAqueous >= system.numPhases()) return 0.0;
            const double amount = state.elementAmountInPhase(data.element, iAqueous);
            const double volume = properties.phaseVolumes()[iAqueous].val;
            const double liter = convertCubicMeterToLiter(volume);
            const double ci = liter ? amount/liter : 0.0;
            return convert_and_apply_scale(ci, data);
        };

        auto species_amount = [=]() -> double
        {
            return convert_and_apply_scale(state.speciesAmount(data.species), data);
        };

        auto species_mass = [=]() -> double
        {
            const double molar_mass = system.species(data.species).molarMass();
            const double amount = state.speciesAmount(data.species);
            return convert_and_apply_scale(amount * molar_mass, data);
        };

        auto species_molality = [=]() -> double
        {
            // Return zero if no water species
            if(iH2O >= system.numSpecies()) return 0.0;
            const double amount = state.speciesAmount(data.species);
            const double kgH2O = state.speciesAmount(iH2O) * waterMolarMass;
            const double mi = kgH2O ? amount/kgH2O : 0.0;
            return convert_and_apply_scale(mi, data);
        };

        auto species_molarity = [=]() -> double
        {
            // Return zero if no aqueous phase
            if(iAqueous >= system.numPhases()) return 0.0;
            const double amount = state.speciesAmount(data.species);
            const double volume = properties.phaseVolumes()[iAqueous].val;
            const double liter = convertCubicMeterToLiter(volume);
            const double ci = liter ? amount/liter : 0.0;
            return convert_and_apply_scale(ci, data);
        };

        auto species_molarfraction = [=]() -> double
        {
            const double xi = properties.molarFractions().val[data.species];
            return convert_and_apply_scale(xi, data);
        };

        auto species_activity = [=]() -> double
        {
            const double ln_ai = properties.lnActivities().val[data.species];
            return convert_and_apply_scale(std::exp(ln_ai), data);
        };

        auto species_activitycoefficient = [=]() -> double
        {
            const double ln_gi = properties.lnActivityCoefficients().val[data.species];
            return convert_and_apply_scale(std::exp(ln_gi), data);
        };

        auto phase_amount = [=]() -> double
        {
            return convert_and_apply_scale(state.phaseAmount(data.phase), data);
        };

        auto phase_mass = [=]() -> double
        {
            const double mass = properties.phaseMasses().val[data.phase];
            return convert_and_apply_scale(mass, data);
        };

        auto ph = [=]() -> double
        {
            // Return zero if no hydron species
            if(iH >= system.numSpecies()) return 0.0;
            const double ln_aH = properties.lnActivities().val[iH];
            return -ln_aH/ln10;
        };

        auto reaction_rate = [=]() -> double
        {
            // Return zero if there are no reactions
            if(r.val.rows() == 0) return 0.0;
            const double ri = r.val[data.reaction];
            return convert_and_apply_scale(ri, data);
        };

        auto reaction_equilibriumindex = [=]() -> double
        {
            // Return zero if there are no reactions
            if(r.val.rows() == 0) return 0.0;
            const double ln_omega = reactions.reaction(data.reaction).lnEquilibriumIndex(properties).val;
            return convert_and_apply_scale(std::exp(ln_omega), data);
        };

        switch(data.quantity)
        {
        case QuantityType::Time: return time;
        case QuantityType::Progress: return progress;
        case QuantityType::Temperature: return temperature;
        case QuantityType::Pressure: return pressure;
        case QuantityType::Amount:
            if(data.element < system.numElements() && data.phase < system.numPhases())
                return element_amount_in_phase;
            if(data.element < system.numElements())
                return element_amount;
            if(data.species < system.numSpecies())
                return species_amount;
            if(data.phase < system.numPhases())
                return phase_amount;
            break;
        case QuantityType::Mass:
            if(data.element < system.numElements() && data.phase < system.numPhases())
                return element_mass_in_phase;
            if(data.element < system.numElements())
                return element_mass;
            if(data.species < system.numSpecies())
                return species_mass;
            if(data.phase < system.numPhases())
                return phase_mass;
            break;
        case QuantityType::Molality:
            if(data.element < system.numElements())
                return element_molality;
            if(data.species < system.numSpecies())
                return species_molality;
            break;
        case QuantityType::Molarity:
            if(data.element < system.numElements())
                return element_molarity;
            if(data.species < system.numSpecies())
                return species_molarity;
            break;
        case QuantityType::Activity: return species_activity;
        case QuantityType::ActivityCoefficient: return species_activitycoefficient;
        case QuantityType::MolarFraction: return species_molarfraction;
        case QuantityType::Fugacity: return species_activity;
        case QuantityType::pH: return ph;
        case QuantityType::Rate: return reaction_rate;
        case QuantityType::EquilibriumIndex: return reaction_equilibriumindex;
        default: break;
        }
        return []() { return 0.0; };
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
