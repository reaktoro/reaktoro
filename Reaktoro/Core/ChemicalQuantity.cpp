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

using ChemicalQuantityFunction = decltype(fn::temperature);

std::map<std::string, ChemicalQuantityFunction> map_quantity_fns =
{
    {"temperature"              , fn::temperature},
    {"pressure"                 , fn::pressure},
    {"volume"                   , fn::volume},
    {"activity"                 , fn::activity},
    {"activitycoefficient"      , fn::activityCoefficient},
    {"fugacity"                 , fn::fugacity},
    {"elementamount"            , fn::elementAmount},
    {"elementamountinphase"     , fn::elementAmountInPhase},
    {"elementmass"              , fn::elementMass},
    {"elementmassinphase"       , fn::elementMassInPhase},
    {"elementmolality"          , fn::elementMolality},
    {"elementmolarity"          , fn::elementMolarity},
    {"speciesamount"            , fn::speciesAmount},
    {"speciesmass"              , fn::speciesMass},
    {"speciesmolefraction"      , fn::speciesMoleFraction},
    {"speciesmolality"          , fn::speciesMolality},
    {"speciesmolarity"          , fn::speciesMolarity},
    {"phaseamount"              , fn::phaseAmount},
    {"phasemass"                , fn::phaseMass},
    {"phasevolume"              , fn::phaseVolume},
    {"ph"                       , fn::pH},
    {"pe"                       , fn::pE},
    {"eh"                       , fn::Eh},
    {"ionicstrength"            , fn::ionicStrength},
    {"fluidvolume"              , fn::fluidVolume},
    {"solidvolume"              , fn::solidVolume},
    {"porosity"                 , fn::porosity},
    {"reactionrate"             , fn::reactionRate},
    {"reactionequilibriumindex" , fn::reactionEquilibriumIndex},
    {"t"                        , fn::t},
    {"time"                     , fn::time},
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

    auto function(std::string str) -> Function
    {
        auto it = function_map.find(str);
        if(it != function_map.end())
            return it->second;

        auto ibegin = str.find_first_of("(");
        auto iend = str.find_last_of(")");

        const std::string funcion_name = str.substr(0, ibegin);
        const std::string arguments = str.substr(ibegin + 1, iend);



        Function newfunc = createFunction(data);

        function_map.insert({str, newfunc});

        return newfunc;
    }

    auto value(std::string str) -> double
    {
        return function(str)();
    }
};

ChemicalQuantity::ChemicalQuantity()
: pimpl(new Impl())
{}

ChemicalQuantity::ChemicalQuantity(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

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

namespace fn {

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
    const ChemicalProperties& properties = quantity.properties();
    const std::string units = args.argument("units", "m3");
    const double factor = units::convert(1.0, units, "m3");
    auto func = [=]() -> double
    {
        const double val = properties.volume().val;
        return factor * val;
    };
    return func;
}

auto activity(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalSystem& system = quantity.system();
    const ChemicalProperties& properties = quantity.properties();
    const std::string species = args.argument(0);
    const Index ispecies = system.indexSpeciesWithError(species);
    auto func = [=]() -> double
    {
        const double ln_ai = properties.lnActivities().val[ispecies];
        return std::exp(ln_ai);
    };
    return func;
}

auto activityCoefficient(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalSystem& system = quantity.system();
    const ChemicalProperties& properties = quantity.properties();
    const std::string species = args.argument(0);
    const Index ispecies = system.indexSpeciesWithError(species);
    auto func = [=]() -> double
    {
        const double ln_gi = properties.lnActivityCoefficients().val[ispecies];
        return std::exp(ln_gi);
    };
    return func;
}

auto fugacity(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalSystem& system = quantity.system();
    const ChemicalProperties& properties = quantity.properties();
    const std::string species = args.argument(0);
    const Index ispecies = system.indexSpeciesWithError(species);
    const std::string units = args.argument("units", "bar");
    const double factor = units::convert(1.0, units, "bar");
    auto func = [=]() -> double
    {
        const double ln_ai = properties.lnActivities().val[ispecies];
        const double val = std::exp(ln_ai);
        return factor * val;
    };
    return func;
}

auto elementAmount(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalState& state = quantity.state();
    const ChemicalSystem& system = quantity.system();
    const std::string element = args.argument(0);
    const std::string units = args.argument("units", "mol");
    const double factor = units::convert(1.0, units, "mol");
    const Index ielement = system.indexElementWithError(element);
    auto func = [=]() -> double
    {
        const double val = state.elementAmount(ielement);
        return factor * val;
    };
    return func;
}

auto elementAmountInPhase(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalState& state = quantity.state();
    const ChemicalSystem& system = quantity.system();
    const std::string element = args.argument(0);
    const std::string phase = args.argument(1);
    const std::string units = args.argument("units", "mol");
    const double factor = units::convert(1.0, units, "mol");
    const Index ielement = system.indexElementWithError(element);
    const Index iphase = system.indexPhaseWithError(phase);
    auto func = [=]() -> double
    {
        const double val = state.elementAmountInPhase(ielement, iphase);
        return factor * val;
    };
    return func;
}

auto elementMass(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalState& state = quantity.state();
    const ChemicalSystem& system = quantity.system();
    const std::string element = args.argument(0);
    const std::string units = args.argument("units", "kg");
    const double factor = units::convert(1.0, units, "kg");
    const Index ielement = system.indexElementWithError(element);
    const double molar_mass = system.element(ielement).molarMass();
    auto func = [=]() -> double
    {
        const double val = state.elementAmount(ielement);
        return factor * molar_mass * val;
    };
    return func;
}

auto elementMassInPhase(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalState& state = quantity.state();
    const ChemicalSystem& system = quantity.system();
    const std::string element = args.argument(0);
    const std::string phase = args.argument(1);
    const std::string units = args.argument("units", "kg");
    const double factor = units::convert(1.0, units, "kg");
    const Index ielement = system.indexElementWithError(element);
    const Index iphase = system.indexPhaseWithError(phase);
    const double molar_mass = system.element(ielement).molarMass();
    auto func = [=]() -> double
    {
        const double val = state.elementAmountInPhase(ielement, iphase);
        return factor * molar_mass * val;
    };
    return func;
}

auto elementMolality(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalState& state = quantity.state();
    const ChemicalSystem& system = quantity.system();
    const std::string element = args.argument(0);
    const std::string units = args.argument("units", "molal");
    const double factor = units::convert(1.0, units, "molal");
    const Index ielement = system.indexElementWithError(element);
    const Index iphase = system.indexPhaseWithError("Aqueous");
    const Index iwater = system.indexSpeciesAny(alternativeWaterNames());
    auto func = [=]() -> double
    {
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
    const ChemicalState& state = quantity.state();
    const ChemicalSystem& system = quantity.system();
    const ChemicalProperties& properties = quantity.properties();
    const std::string element = args.argument(0);
    const std::string units = args.argument("units", "molar");
    const double factor = units::convert(1.0, units, "molar");
    const Index ielement = system.indexElementWithError(element);
    const Index iphase = system.indexPhaseWithError("Aqueous");
    auto func = [=]() -> double
    {
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
    const ChemicalState& state = quantity.state();
    const ChemicalSystem& system = quantity.system();
    const std::string species = args.argument(0);
    const Index ispecies = system.indexSpeciesWithError(species);
    const std::string units = args.argument("units", "mol");
    const double factor = units::convert(1.0, units, "mol");
    auto func = [=]() -> double
    {
        const double val = state.speciesAmount(ispecies);
        return factor * val;
    };
    return func;
}

auto speciesMass(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalState& state = quantity.state();
    const ChemicalSystem& system = quantity.system();
    const std::string species = args.argument(0);
    const Index ispecies = system.indexSpeciesWithError(species);
    const double molar_mass = system.species(ispecies).molarMass();
    const std::string units = args.argument("units", "kg");
    const double factor = units::convert(1.0, units, "kg");
    auto func = [=]() -> double
    {
        const double val = state.speciesAmount(ispecies);
        return factor * molar_mass * val;
    };
    return func;
}

auto speciesMoleFraction(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalProperties& properties = quantity.properties();
    const ChemicalSystem& system = quantity.system();
    const std::string species = args.argument(0);
    const Index ispecies = system.indexSpeciesWithError(species);
    auto func = [=]() -> double
    {
        const double xi = properties.molarFractions().val[ispecies];
        return xi;
    };
    return func;
}

auto speciesMolality(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalState& state = quantity.state();
    const ChemicalSystem& system = quantity.system();
    const std::string species = args.argument(0);
    const std::string units = args.argument("units", "molal");
    const double factor = units::convert(1.0, units, "molal");
    const Index ispecies = system.indexSpeciesWithError(species);
    const Index iwater = system.indexSpeciesAny(alternativeWaterNames());
    auto func = [=]() -> double
    {
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
    const ChemicalState& state = quantity.state();
    const ChemicalSystem& system = quantity.system();
    const ChemicalProperties& properties = quantity.properties();
    const std::string species = args.argument(0);
    const std::string units = args.argument("units", "molar");
    const double factor = units::convert(1.0, units, "molar");
    const Index ispecies = system.indexSpeciesWithError(species);
    const Index iphase = system.indexPhaseWithError("Aqueous");
    auto func = [=]() -> double
    {
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
    const ChemicalProperties& properties = quantity.properties();
    const ChemicalSystem& system = quantity.system();
    const std::string phase = args.argument(0);
    const Index iphase = system.indexPhaseWithError(phase);
    const std::string units = args.argument("units", "mol");
    const double factor = units::convert(1.0, units, "mol");
    auto func = [=]() -> double
    {
        const double val = properties.phaseAmounts().val[iphase];
        return factor * val;
    };
    return func;
}

auto phaseMass(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalProperties& properties = quantity.properties();
    const ChemicalSystem& system = quantity.system();
    const std::string phase = args.argument(0);
    const Index iphase = system.indexPhaseWithError(phase);
    const std::string units = args.argument("units", "kg");
    const double factor = units::convert(1.0, units, "kg");
    auto func = [=]() -> double
    {
        const double val = properties.phaseMasses().val[iphase];
        return factor * val;
    };
    return func;
}

auto phaseVolume(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalProperties& properties = quantity.properties();
    const ChemicalSystem& system = quantity.system();
    const std::string phase = args.argument(0);
    const Index iphase = system.indexPhaseWithError(phase);
    const std::string units = args.argument("units", "m3");
    const double factor = units::convert(1.0, units, "m3");
    auto func = [=]() -> double
    {
        const double val = properties.phaseVolumes().val[iphase];
        return factor * val;
    };
    return func;
}

auto pH(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalProperties& properties = quantity.properties();
    auto func = [=]() -> double
    {
        const double val = properties.aqueous().pH().val;
        return val;
    };
    return func;
}

auto pE(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalProperties& properties = quantity.properties();
    auto func = [=]() -> double
    {
        const double val = properties.aqueous().pE().val;
        return val;
    };
    return func;
}

auto Eh(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalProperties& properties = quantity.properties();
    const std::string units = args.argument("units", "volt");
    const double factor = units::convert(1.0, units, "volt");
    auto func = [=]() -> double
    {
        const double val = properties.aqueous().Eh().val;
        return factor * val;
    };
    return func;
}

auto ionicStrength(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalProperties& properties = quantity.properties();
    const std::string units = args.argument("units", "molal");
    const double factor = units::convert(1.0, units, "molal");
    auto func = [=]() -> double
    {
        const double val = properties.aqueous().ionicStrength().val;
        return factor * val;
    };
    return func;
}

auto fluidVolume(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalProperties& properties = quantity.properties();
    const std::string units = args.argument("units", "m3");
    const double factor = units::convert(1.0, units, "m3");
    auto func = [=]() -> double
    {
        const double val = properties.fluidVolume().val;
        return factor * val;
    };
    return func;
}

auto solidVolume(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ChemicalProperties& properties = quantity.properties();
    const std::string units = args.argument("units", "m3");
    const double factor = units::convert(1.0, units, "m3");
    auto func = [=]() -> double
    {
        const double val = properties.solidVolume().val;
        return factor * val;
    };
    return func;
}

auto reactionRate(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ReactionSystem& reactions = quantity.reactions();
    const ChemicalVector& rates = quantity.rates();
    const std::string reaction = args.argument(0);
    const Index ireaction = reactions.indexReactionWithError(reaction);
    const std::string units = args.argument("units", "mol/s");
    const double factor = units::convert(1.0, units, "mol/s");
    auto func = [=]() -> double
    {
        const double val = rates.val[ireaction];
        return factor * val;
    };
    return func;
}

auto reactionEquilibriumIndex(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const ReactionSystem& reactions = quantity.reactions();
    const ChemicalProperties& properties = quantity.properties();
    const std::string reaction = args.argument(0);
    const Index ireaction = reactions.indexReactionWithError(reaction);
    auto func = [=]() -> double
    {
        const double ln_omega = reactions.reaction(ireaction).lnEquilibriumIndex(properties).val;
        return std::exp(ln_omega);
    };
    return func;
}

auto t(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    return time(quantity, arguments);
}

auto time(const ChemicalQuantity& quantity, std::string arguments) -> std::function<double()>
{
    const Args args(arguments);
    const std::string units = args.argument("units", "s");
    const double factor = units::convert(1.0, units, "s");
    auto func = [=]() -> double
    {
        const double val = quantity.tag();
        return factor * val;
    };
    return func;
}

} // namespace fn
} // namespace Reaktoro
