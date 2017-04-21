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
namespace quantity {

auto temperature(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto pressure(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto volume(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto moleFraction(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto activity(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto activityCoefficient(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto fugacity(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto chemicalPotential(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto elementAmount(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto elementAmountInPhase(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto elementMass(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto elementMassInPhase(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto elementMolality(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto elementMolarity(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto speciesAmount(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto speciesMass(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto speciesMolality(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto speciesMolarity(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto phaseAmount(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto phaseMass(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto phaseVolume(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto pH(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto pE(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto Eh(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto ionicStrength(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto fluidVolume(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto solidVolume(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto porosity(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto reactionRate(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto reactionEquilibriumIndex(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

auto tag(const ChemicalQuantity& quantity, std::string args) -> std::function<double()>;

using Function = std::function<decltype(temperature)>;

const std::map<std::string, Function> fndict =
{
    {"temperature"              , quantity::temperature},
    {"pressure"                 , quantity::pressure},
    {"volume"                   , quantity::volume},
    {"molefraction"             , quantity::moleFraction},
    {"activity"                 , quantity::activity},
    {"activitycoefficient"      , quantity::activityCoefficient},
    {"fugacity"                 , quantity::fugacity},
    {"chemicalpotential"        , quantity::chemicalPotential},
    {"elementamount"            , quantity::elementAmount},
    {"elementamountinphase"     , quantity::elementAmountInPhase},
    {"elementmass"              , quantity::elementMass},
    {"elementmassinphase"       , quantity::elementMassInPhase},
    {"elementmolality"          , quantity::elementMolality},
    {"elementmolarity"          , quantity::elementMolarity},
    {"speciesamount"            , quantity::speciesAmount},
    {"speciesmass"              , quantity::speciesMass},
    {"speciesmolality"          , quantity::speciesMolality},
    {"speciesmolarity"          , quantity::speciesMolarity},
    {"phaseamount"              , quantity::phaseAmount},
    {"phasemass"                , quantity::phaseMass},
    {"phasevolume"              , quantity::phaseVolume},
    {"ph"                       , quantity::pH},
    {"pe"                       , quantity::pE},
    {"eh"                       , quantity::Eh},
    {"ionicstrength"            , quantity::ionicStrength},
    {"fluidvolume"              , quantity::fluidVolume},
    {"solidvolume"              , quantity::solidVolume},
    {"porosity"                 , quantity::porosity},
    {"reactionrate"             , quantity::reactionRate},
    {"reactionequilibriumindex" , quantity::reactionEquilibriumIndex},
    {"tag"                      , quantity::tag},
    {"t"                        , quantity::tag},
    {"time"                     , quantity::tag},
    {"progress"                 , quantity::tag},
};

auto function(std::string fname) -> Function;

} // namespace quantity

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
    double tag;

    /// The temperature of the chemical system (in units of K).
    double T;

    /// The pressure of the chemical system (in units of Pa).
    double P;

    /// The molar amounts of the species in the chemical system (in units of mol).
    Vector n;

    /// The rates of the reactions in the chemical system (in units of mol/s).
    ChemicalVector rates;

    /// All created chemical quantity functions from formatted strings
    std::map<std::string, Function> function_map;

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a custom Impl instance with given ChemicalSystem object
    Impl(const ChemicalSystem& system)
    : system(system)
    {
    }

    /// Construct a custom Impl instance with given ReactionSystem object
    Impl(const ReactionSystem& reactions)
    : system(reactions.system()), reactions(reactions)
    {
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
        tag = t_;

        // Update the temperature, pressure and molar composition of the system
        T = state.temperature();
        P = state.pressure();
        n = state.speciesAmounts();

        // Update the thermodynamic properties of the system
        properties = system.properties(T, P, n);

        // Update the rates of the reactions
        if(!reactions.reactions().empty())
            rates = reactions.rates(properties);
    }

    auto function(const ChemicalQuantity& quantity, std::string str) -> Function
    {
        auto it = function_map.find(str);
        if(it != function_map.end())
            return it->second;

        str = trim(str);
        auto ibegin = str.find_first_of("(");
        auto iend = str.find_last_of(")");

        const std::string fname = lowercase(str.substr(0, ibegin));
        const std::string arguments = ibegin == std::string::npos ?
                "" : str.substr(ibegin+1, iend-ibegin-1);

        Function newfunc = quantity::function(fname)(quantity, arguments);

        function_map.insert({str, newfunc});

        return newfunc;
    }

    auto value(const ChemicalQuantity& quantity, std::string str) -> double
    {
        return function(quantity, str)();
    }
};

ChemicalQuantity::ChemicalQuantity()
: pimpl(new Impl())
{}

ChemicalQuantity::ChemicalQuantity(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

ChemicalQuantity::ChemicalQuantity(const ChemicalState& state)
: ChemicalQuantity(state.system())
{
    update(state);
}

ChemicalQuantity::ChemicalQuantity(const ReactionSystem& reactions)
: pimpl(new Impl(reactions))
{}

ChemicalQuantity::~ChemicalQuantity()
{}

auto ChemicalQuantity::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto ChemicalQuantity::reactions() const -> const ReactionSystem&
{
    return pimpl->reactions;
}

auto ChemicalQuantity::state() const -> const ChemicalState&
{
    return pimpl->state;
}

auto ChemicalQuantity::properties() const -> const ChemicalProperties&
{
    return pimpl->properties;
}

auto ChemicalQuantity::rates() const -> const ChemicalVector&
{
    return pimpl->rates;
}

auto ChemicalQuantity::tag() const -> double
{
    return pimpl->tag;
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
    return pimpl->value(*this, str);
}

auto ChemicalQuantity::function(std::string str) const -> Function
{
    return pimpl->function(*this, str);
}

auto ChemicalQuantity::operator[](std::string quantity) const -> double
{
    return value(quantity);
}

namespace quantity {

/// A type used to describe the list of arguments for quantity querying.
struct Args
{
    /// Construct a default Args instance.
    Args()
    {}

    /// Construct a custom Args instance from a formatted string.
    Args(std::string arguments)
    : str(arguments)
    {
        auto words = split(str);
        for(auto word : words)
        {
            auto pair = split(word, "=");
            if(pair.size() == 2)
                kwargs.insert({ pair[0], pair[1] });
            else args.push_back(word);
        }
    }

    /// Return the ith non-keyword argument.
    auto argument(Index i) const -> std::string
    {
        assert(i < args.size());
        return args[i];
    }

    /// Return the argument with given keyword.
    auto argument(std::string keyword, std::string ifnone = "") const -> std::string
    {
        auto iter = kwargs.find(keyword);
        if(iter != kwargs.end()) return iter->second;
        else return ifnone;
    }

private:
    /// The original formated string
    std::string str;

    /// The non-keyword arguments
    std::vector<std::string> args;

    /// The keyword arguments
    std::map<std::string, std::string> kwargs;
};

auto temperature(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const std::string units = args.argument("units", "K");
    auto func = [=]() -> double
    {
        const double val = quantity.state().temperature();
        return units::convert(val, "K", units);
    };
    return func;
}

auto pressure(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const std::string units = args.argument("units", "Pa");
    auto func = [=]() -> double
    {
        const double val = quantity.state().pressure();
        return units::convert(val, "Pa", units);
    };
    return func;
}

auto volume(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const std::string units = args.argument("units", "m3");
    const double factor = units::convert(1.0, "m3", units);
    auto func = [=]() -> double
    {
        const ChemicalProperties& properties = quantity.properties();
        const double val = properties.volume().val;
        return factor * val;
    };
    return func;
}

auto moleFraction(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalSystem& system = quantity.system();
    const std::string species = args.argument(0);
    const Index ispecies = system.indexSpeciesWithError(species);
    auto func = [=]() -> double
    {
        const ChemicalProperties& properties = quantity.properties();
        const double xi = properties.molarFractions().val[ispecies];
        return xi;
    };
    return func;
}

auto activity(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalSystem& system = quantity.system();
    const std::string species = args.argument(0);
    const Index ispecies = system.indexSpeciesWithError(species);
    auto func = [=]() -> double
    {
        const ChemicalProperties& properties = quantity.properties();
        const double ln_ai = properties.lnActivities().val[ispecies];
        return std::exp(ln_ai);
    };
    return func;
}

auto activityCoefficient(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalSystem& system = quantity.system();
    const std::string species = args.argument(0);
    const Index ispecies = system.indexSpeciesWithError(species);
    auto func = [=]() -> double
    {
        const ChemicalProperties& properties = quantity.properties();
        const double ln_gi = properties.lnActivityCoefficients().val[ispecies];
        return std::exp(ln_gi);
    };
    return func;
}

auto fugacity(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalSystem& system = quantity.system();
    const std::string species = args.argument(0);
    const Index ispecies = system.indexSpeciesWithError(species);
    const std::string units = args.argument("units", "bar");
    const double factor = units::convert(1.0, "bar", units);
    auto func = [=]() -> double
    {
        const ChemicalProperties& properties = quantity.properties();
        const double ln_ai = properties.lnActivities().val[ispecies];
        const double val = std::exp(ln_ai);
        return factor * val;
    };
    return func;
}

auto chemicalPotential(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalSystem& system = quantity.system();
    const std::string species = args.argument(0);
    const Index ispecies = system.indexSpeciesWithError(species);
    const std::string units = args.argument("units", "J/mol");
    const double factor = units::convert(1.0, "J/mol", units);
    auto func = [=]() -> double
    {
        const ChemicalProperties& properties = quantity.properties();
        const double val = properties.chemicalPotentials().val[ispecies];
        return factor * val;
    };
    return func;
}

auto elementAmount(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalSystem& system = quantity.system();
    const std::string element = args.argument(0);
    const std::string units = args.argument("units", "mol");
    const double factor = units::convert(1.0, "mol", units);
    const Index ielement = system.indexElementWithError(element);
    auto func = [=]() -> double
    {
        const ChemicalState& state = quantity.state();
        const double val = state.elementAmount(ielement);
        return factor * val;
    };
    return func;
}

auto elementAmountInPhase(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalSystem& system = quantity.system();
    const std::string element = args.argument(0);
    const std::string phase = args.argument(1);
    const std::string units = args.argument("units", "mol");
    const double factor = units::convert(1.0, "mol", units);
    const Index ielement = system.indexElementWithError(element);
    const Index iphase = system.indexPhaseWithError(phase);
    auto func = [=]() -> double
    {
        const ChemicalState& state = quantity.state();
        const double val = state.elementAmountInPhase(ielement, iphase);
        return factor * val;
    };
    return func;
}

auto elementMass(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalSystem& system = quantity.system();
    const std::string element = args.argument(0);
    const std::string units = args.argument("units", "kg");
    const double factor = units::convert(1.0, "kg", units);
    const Index ielement = system.indexElementWithError(element);
    const double molar_mass = system.element(ielement).molarMass();
    auto func = [=]() -> double
    {
        const ChemicalState& state = quantity.state();
        const double val = state.elementAmount(ielement);
        return factor * molar_mass * val;
    };
    return func;
}

auto elementMassInPhase(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalSystem& system = quantity.system();
    const std::string element = args.argument(0);
    const std::string phase = args.argument(1);
    const std::string units = args.argument("units", "kg");
    const double factor = units::convert(1.0, "kg", units);
    const Index ielement = system.indexElementWithError(element);
    const Index iphase = system.indexPhaseWithError(phase);
    const double molar_mass = system.element(ielement).molarMass();
    auto func = [=]() -> double
    {
        const ChemicalState& state = quantity.state();
        const double val = state.elementAmountInPhase(ielement, iphase);
        return factor * molar_mass * val;
    };
    return func;
}

auto elementMolality(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalSystem& system = quantity.system();
    const std::string element = args.argument(0);
    const std::string units = args.argument("units", "molal");
    const double factor = units::convert(1.0, "molal", units);
    const Index ielement = system.indexElementWithError(element);
    const Index iphase = system.indexPhaseWithError("Aqueous");
    const Index iwater = system.indexSpeciesAny(alternativeWaterNames());
    auto func = [=]() -> double
    {
        const ChemicalState& state = quantity.state();
        const double amount = state.elementAmountInPhase(ielement, iphase);
        const double kgH2O = state.speciesAmount(iwater) * waterMolarMass;
        const double mi = kgH2O ? amount/kgH2O : 0.0;
        return factor * mi;
    };
    return func;
}

auto elementMolarity(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalSystem& system = quantity.system();
    const std::string element = args.argument(0);
    const std::string units = args.argument("units", "molar");
    const double factor = units::convert(1.0, "molar", units);
    const Index ielement = system.indexElementWithError(element);
    const Index iphase = system.indexPhaseWithError("Aqueous");
    auto func = [=]() -> double
    {
        const ChemicalProperties& properties = quantity.properties();
        const ChemicalState& state = quantity.state();
        const double amount = state.elementAmountInPhase(ielement, iphase);
        const double volume = properties.phaseVolumes()[iphase].val;
        const double liter = convertCubicMeterToLiter(volume);
        const double ci = liter ? amount/liter : 0.0;
        return factor * ci;
    };
    return func;
}

auto speciesAmount(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalSystem& system = quantity.system();
    const std::string species = args.argument(0);
    const Index ispecies = system.indexSpeciesWithError(species);
    const std::string units = args.argument("units", "mol");
    const double factor = units::convert(1.0, "mol", units);
    auto func = [=]() -> double
    {
        const ChemicalState& state = quantity.state();
        const double val = state.speciesAmount(ispecies);
        return factor * val;
    };
    return func;
}

auto speciesMass(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalSystem& system = quantity.system();
    const std::string species = args.argument(0);
    const Index ispecies = system.indexSpeciesWithError(species);
    const double molar_mass = system.species(ispecies).molarMass();
    const std::string units = args.argument("units", "kg");
    const double factor = units::convert(1.0, "kg", units);
    auto func = [=]() -> double
    {
        const ChemicalState& state = quantity.state();
        const double val = state.speciesAmount(ispecies);
        return factor * molar_mass * val;
    };
    return func;
}

auto speciesMolality(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalSystem& system = quantity.system();
    const std::string species = args.argument(0);
    const std::string units = args.argument("units", "molal");
    const double factor = units::convert(1.0, "molal", units);
    const Index ispecies = system.indexSpeciesWithError(species);
    const Index iwater = system.indexSpeciesAny(alternativeWaterNames());
    auto func = [=]() -> double
    {
        const ChemicalState& state = quantity.state();
        const double amount = state.speciesAmount(ispecies);
        const double kgH2O = state.speciesAmount(iwater) * waterMolarMass;
        const double mi = kgH2O ? amount/kgH2O : 0.0;
        return factor * mi;
    };
    return func;
}

auto speciesMolarity(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalSystem& system = quantity.system();
    const std::string species = args.argument(0);
    const std::string units = args.argument("units", "molar");
    const double factor = units::convert(1.0, "molar", units);
    const Index ispecies = system.indexSpeciesWithError(species);
    const Index iphase = system.indexPhaseWithError("Aqueous");
    auto func = [=]() -> double
    {
        const ChemicalProperties& properties = quantity.properties();
        const ChemicalState& state = quantity.state();
        const double amount = state.speciesAmount(ispecies);
        const double volume = properties.phaseVolumes()[iphase].val;
        const double liter = convertCubicMeterToLiter(volume);
        const double ci = liter ? amount/liter : 0.0;
        return factor * ci;
    };
    return func;
}

auto phaseAmount(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalSystem& system = quantity.system();
    const std::string phase = args.argument(0);
    const Index iphase = system.indexPhaseWithError(phase);
    const std::string units = args.argument("units", "mol");
    const double factor = units::convert(1.0, "mol", units);
    auto func = [=]() -> double
    {
        const ChemicalProperties& properties = quantity.properties();
        const double val = properties.phaseAmounts().val[iphase];
        return factor * val;
    };
    return func;
}

auto phaseMass(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalSystem& system = quantity.system();
    const std::string phase = args.argument(0);
    const Index iphase = system.indexPhaseWithError(phase);
    const std::string units = args.argument("units", "kg");
    const double factor = units::convert(1.0, "kg", units);
    auto func = [=]() -> double
    {
        const ChemicalProperties& properties = quantity.properties();
        const double val = properties.phaseMasses().val[iphase];
        return factor * val;
    };
    return func;
}

auto phaseVolume(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalSystem& system = quantity.system();
    const std::string phase = args.argument(0);
    const Index iphase = system.indexPhaseWithError(phase);
    const std::string units = args.argument("units", "m3");
    const double factor = units::convert(1.0, "m3", units);
    auto func = [=]() -> double
    {
        const ChemicalProperties& properties = quantity.properties();
        const double val = properties.phaseVolumes().val[iphase];
        return factor * val;
    };
    return func;
}

auto pH(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    auto func = [=]() -> double
    {
        const ChemicalProperties& properties = quantity.properties();
        const double val = properties.aqueous().pH().val;
        return val;
    };
    return func;
}

auto pE(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    auto func = [=]() -> double
    {
        const ChemicalProperties& properties = quantity.properties();
        const double val = properties.aqueous().pE().val;
        return val;
    };
    return func;
}

auto Eh(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const std::string units = args.argument("units", "volt");
    const double factor = units::convert(1.0, "volt", units);
    auto func = [=]() -> double
    {
        const ChemicalProperties& properties = quantity.properties();
        const double val = properties.aqueous().Eh().val;
        return factor * val;
    };
    return func;
}

auto ionicStrength(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const std::string units = args.argument("units", "molal");
    const double factor = units::convert(1.0, "molal", units);
    auto func = [=]() -> double
    {
        const ChemicalProperties& properties = quantity.properties();
        const double val = properties.aqueous().ionicStrength().val;
        return factor * val;
    };
    return func;
}

auto fluidVolume(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const std::string units = args.argument("units", "m3");
    const double factor = units::convert(1.0, "m3", units);
    auto func = [=]() -> double
    {
        const ChemicalProperties& properties = quantity.properties();
        const double val = properties.fluidVolume().val;
        return factor * val;
    };
    return func;
}

auto solidVolume(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const std::string units = args.argument("units", "m3");
    const double factor = units::convert(1.0, "m3", units);
    auto func = [=]() -> double
    {
        const ChemicalProperties& properties = quantity.properties();
        const double val = properties.solidVolume().val;
        return factor * val;
    };
    return func;
}

auto porosity(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    auto func = [=]() -> double
    {
        const ChemicalProperties& properties = quantity.properties();
        const double volume = properties.volume().val;
        const double solid_volume = properties.solidVolume().val;
        return 1.0 - solid_volume/volume;
    };
    return func;
}

auto reactionRate(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ReactionSystem& reactions = quantity.reactions();
    const std::string reaction = args.argument(0);
    const Index ireaction = reactions.indexReactionWithError(reaction);
    const std::string units = args.argument("units", "mol/s");
    const double factor = units::convert(1.0, "mol/s", units);
    auto func = [=]() -> double
    {
        const ChemicalVector& rates = quantity.rates();
        const double val = rates.val[ireaction];
        return factor * val;
    };
    return func;
}

auto reactionEquilibriumIndex(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ReactionSystem& reactions = quantity.reactions();
    const std::string reaction = args.argument(0);
    const Index ireaction = reactions.indexReactionWithError(reaction);
    auto func = [=]() -> double
    {
        const ChemicalProperties& properties = quantity.properties();
        const double ln_omega = reactions.reaction(ireaction).lnEquilibriumIndex(properties).val;
        return std::exp(ln_omega);
    };
    return func;
}

auto tag(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const std::string units = args.argument("units", "s");
    const double factor = units::convert(1.0, "s", units);
    auto func = [=]() -> double
    {
        const double val = quantity.tag();
        return factor * val;
    };
    return func;
}

auto function(std::string fname) -> Function
{
    auto iter = quantity::fndict.find(fname);
    Assert(iter != quantity::fndict.end(),
        "Could not create the quantity function with name `" + fname + "`.",
        "This function name has been misspelled or it is not supported.");
    return iter->second;
}

} // namespace quantity
} // namespace Reaktoro
