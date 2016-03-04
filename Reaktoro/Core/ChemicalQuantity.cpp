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

/// The stored result of ln(10)
const double ln10 = 2.30258509299;

/// A type used to describe the chemical quantity type.
enum QuantityType
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
enum OperatorType
{
    None, ln, log, exp
};

/// A type used to describe a chemical quantity.
struct Description
{
    /// The type of the quantity.
    QuantityType quantity = NotSupported;

    /// The index of the element for which the quantity is related.
    /// This can be empty if the quantity is not related to an element.
    Index element;

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

    /// The progress variable at which the chemical state is referred (in units of s)
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
        // Check the string is not empty
        Assert(!str.empty(), "Could not parse the given quantity string.",
            "The quantity string must be non-empty.");

        // Auxiliary strings that are components of the formatted string `str`
        std::string op;
        std::string units;
        std::string entity;
        std::string phase;

        // Extract the operator, if any, and set str = x, where x comes from `op(x)`
        if(str.substr(0, 2) == "ln")  { op = "ln" ; str = str.substr(3, str.size()-4); } // ln(molal(CO2))
        if(str.substr(0, 3) == "log") { op = "log"; str = str.substr(4, str.size()-5); }
        if(str.substr(0, 3) == "exp") { op = "exp"; str = str.substr(4, str.size()-5); }

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
            {"xi"                  , QuantityType::Progress},
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
        if(operators.count(op))
            desc.op = operators.at(op);
        else RuntimeError("Could not identify operator type from `" + op + "`.",
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
