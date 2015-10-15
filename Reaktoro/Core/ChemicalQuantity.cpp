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

#include "ChemicalQuantity.hpp"

// C++ includes
#include <map>
#include <set>

// Reaktoro includes
#include <Reaktoro/Common/ConvertUtils.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {
namespace {

auto defaultQuantityUnits(std::string quantity) -> std::string
{
	static const std::map<std::string, std::string> default_units =
	{
		{"t",        "s"},
		{"time",     "s"},
		{"amount",   "mol"},
		{"mass",     "kg"},
		{"molality", "molal"},
		{"molarity", "molar"},
		{"rate",     "mol/s"},
	};
	auto iter = default_units.find(quantity);
	return iter != default_units.end() ? iter->second : "";
}

struct QuantityInfo
{
	/// The name of the quantity (e.g., `amount`, `mass`, `molality`, `activity`, etc.)
	std::string quantity;

	/// The name of the element for which the quantity is related.
	/// This can be empty if the quantity is not related to an element.
	/// For example, `mass[Quartz]:mol` results in `quantity = "mass"`,
	/// `species = "Quartz"`, and `units = "mol". The member `element` is then empty.
	std::string element;

	/// The name of the species for which the quantity is related.
	/// This can be empty if the quantity is not related to a species.
	/// For example, `amount[H]:mmol` results in `quantity = "amount"`,
	/// `element = "H"`, and `units = "mmol". The member `species` is then empty.
	std::string species;

	/// The name of the phase for which the quantity is related.
	/// This can be empty if the quantity is not related to a phase.
	/// An example in which it is not empty is `amount[Gaseous]`, which results in
	/// `quantity = "amount"`, and `phase = "Gaseous"`. In this case, members
	/// `element` and `species` are empty.
	std::string phase;

	/// The name of the reaction for which the quantity is related.
	/// Reaction related quantities include `rate` and `equilibriumIndex`,
	/// and `reactionQuotient`.
	std::string reaction;

	/// The name of the units for the quantity
	std::string units;

	/// The name of the scale for the quantity (log, ln, exp)
	std::string scale;
};

auto parseQuantityString(std::string str) -> QuantityInfo
{
	Assert(!str.empty(), "Could not parse the quantity string `" + str + "`.",
		"The quantity string must be non-empty.");

	QuantityInfo info;

	// Split the string into individual words at the spaces
	auto words = split(str);

	// Set the name of the quantity (the first word in a list of space-separated words)
	info.quantity = words[0];

	// Iterate over keyword-value pairs (e.g., element=H, species=CO2(aq), units=mol/s)
	for(unsigned i = 1; i < words.size(); ++i)
	{
		auto pair = split(words[i], "=");

		Assert(pair.size() == 2, "Could not parse the quantity string `" + str + "`.",
			"The keyword-value pair `" + words[i] + "` is invalid. "
				"Do not use space between `=`, e.g., use species=H2O(l), not species = H2O(l).");

		if(pair[0] == "element")
		{
			Assert(info.element.empty(), "Could not parse the quantity string `" + str + "`.",
				"The keyword `element` appears more than once in the quantity string.");
			info.element = pair[1];
		}
		else if(pair[0] == "species")
		{
			Assert(info.species.empty(), "Could not parse the quantity string `" + str + "`.",
				"The keyword `species` appears more than once in the quantity string.");
			info.species = pair[1];
		}
		else if(pair[0] == "phase")
		{
			Assert(info.phase.empty(), "Could not parse the quantity string `" + str + "`.",
				"The keyword `phase` appears more than once in the quantity string.");
			info.phase = pair[1];
		}
		else if(pair[0] == "reaction")
		{
			Assert(info.reaction.empty(), "Could not parse the quantity string `" + str + "`.",
				"The keyword `reaction` appears more than once in the quantity string.");
			info.reaction = pair[1];
		}
		else if(pair[0] == "units")
		{
			Assert(info.units.empty(), "Could not parse the quantity string `" + str + "`.",
				"The keyword `units` appears more than once in the quantity string.");
			info.units = pair[1];
		}
		else if(pair[0] == "scale")
		{
			Assert(info.scale.empty(), "Could not parse the quantity string `" + str + "`.",
				"The keyword `scale` appears more than once in the quantity string.");
			info.scale = pair[1];
		}
	}

	// Ensure the default quantity units is set in case none was specified
	info.units = info.units.empty() ? defaultQuantityUnits(info.quantity) : info.units;

	return info;
}

auto validateQuantityInfo(const QuantityInfo& info, const ChemicalSystem& system, const ReactionSystem& reactions) -> void
{
	/// The valid quantity names
	static const std::set<std::string> valid_quantity_names =
	{
		"t", "time", "amount", "mass", "molality", "molarity", "activity", "pH",
		"activityCoefficient", "molarFraction", "rate", "equilibriumIndex"
	};

	/// The valid scale names
	static const std::set<std::string> valid_scale_names =
	{
		"log", "ln", "exp"
	};

	/// The quantities that accept element keyword
	static const std::set<std::string> quantities_element_ok =
	{
		"amount", "mass", "molality", "molarity"
	};

	/// The quantities that accept species keyword
	static const std::set<std::string> quantities_species_ok =
	{
		"amount", "mass", "molality", "molarity", "activity", "activityCoefficient", "molarFraction"
	};

	/// The quantities that accept phase keyword
	static const std::set<std::string> quantities_phase_ok =
	{
		"amount", "mass"
	};

	/// The quantities that accept reaction keyword
	static const std::set<std::string> quantities_reaction_ok =
	{
		"rate", "equilibriumIndex"
	};

	/// The quantities that accept units keyword
	static const std::set<std::string> quantities_units_ok =
	{
		"t", "time", "amount", "mass", "molality", "molarity", "rate"
	};

	/// The quantities that require species keyword
	static const std::set<std::string> quantities_require_species =
	{
		"activity", "activityCoefficient", "molarFraction"
	};

	/// The quantities that require reaction keyword
	static const std::set<std::string> quantities_require_reaction =
	{
		"rate", "equilibriumIndex"
	};

	Assert(valid_quantity_names.count(info.quantity),
		"Could not validate the quantity string.",
			"The quantity name `" + info.quantity + "` is not supported.");

	if(info.element.size())
		Assert(system.indexElement(info.element) < system.numElements(),
			"Could not validate the quantity string.",
				"The specified element `" + info.element + "` is not present in the chemical system.");

	if(info.species.size())
		Assert(system.indexSpecies(info.species) < system.numSpecies(),
			"Could not validate the quantity string.",
				"The specified species `" + info.species + "` is not present in the chemical system.");

	if(info.phase.size())
		Assert(system.indexPhase(info.phase) < system.numPhases(),
			"Could not validate the quantity string.",
				"The specified phase `" + info.phase + "` is not present in the chemical system.");

	if(info.reaction.size())
		Assert(reactions.indexReaction(info.reaction) < reactions.numReactions(),
			"Could not validate the quantity string.",
				"The specified reaction `" + info.reaction + "` is not present in the reaction system.");

	if(info.scale.size())
		Assert(valid_scale_names.count(info.scale),
			"Could not validate the quantity string.",
				"The specified scale `" + info.scale + "` is not supported.");

	if(info.element.size())
		Assert(quantities_element_ok.count(info.quantity),
			"Could not validate the quantity string.",
				"The specified quantity `" + info.quantity + "` does not accept `element` keyword.");

	if(info.species.size())
		Assert(quantities_species_ok.count(info.quantity),
			"Could not validate the quantity string.",
				"The specified quantity `" + info.quantity + "` does not accept `species` keyword.");

	if(info.phase.size())
		Assert(quantities_phase_ok.count(info.quantity),
			"Could not validate the quantity string.",
				"The specified quantity `" + info.quantity + "` does not accept `phase` keyword.");

	if(info.reaction.size())
		Assert(quantities_reaction_ok.count(info.quantity),
			"Could not validate the quantity string.",
				"The specified quantity `" + info.quantity + "` does not accept `reaction` keyword.");

	if(info.units.size())
		Assert(quantities_units_ok.count(info.quantity),
			"Could not validate the quantity string.",
				"The specified quantity `" + info.quantity + "` does not accept `units` keyword.");

	if(quantities_require_species.count(info.quantity))
		Assert(info.species.size(),
			"Could not validate the quantity string.",
				"The specified quantity `" + info.quantity + "` requires a `species` keyword.");

	if(quantities_require_reaction.count(info.quantity))
		Assert(info.reaction.size(),
			"Could not validate the quantity string.",
				"The specified quantity `" + info.quantity + "` requires a `reaction` keyword.");
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

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a custom Impl instance with given ChemicalSystem object
    Impl(const ChemicalSystem& system)
    : system(system)
    {}

    /// Construct a custom Impl instance with given ReactionSystem object
    Impl(const ReactionSystem& reactions)
    : system(reactions.system()), reactions(reactions)
    {}

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

    auto value(std::string str) const -> double
    {
    	QuantityInfo info = parseQuantityString(str);

    	validateQuantityInfo(info, system, reactions);

        if(info.quantity == "t" || info.quantity == "time")
        {
            return units::convert(t, "s", info.units);
        }
        if(info.quantity == "amount")
        {
        	if(!info.species.empty())
				return state.speciesAmount(info.species, info.units);
        	if(!info.element.empty() && info.phase.empty())
				return state.elementAmount(info.element, info.units);
        	if(!info.element.empty() && !info.phase.empty())
				return state.elementAmountInPhase(info.element, info.phase, info.units);
        	if(!info.phase.empty() && info.element.empty())
				return state.phaseAmount(info.phase, info.units);
        }
        if(info.quantity == "molarFraction")
        {
            auto ispecies = system.indexSpeciesWithError(info.species);
            return properties.molarFractions()[ispecies].val;
        }
        if(info.quantity == "molality")
        {
        	double amount = 0.0;
        	if(info.species.size())
				amount = state.speciesAmount(info.species);
        	if(info.element.size())
				amount = state.elementAmountInPhase(info.element, "Aqueous");
            const double kgH2O = state.speciesAmount("H2O(l)") * waterMolarMass;
            const double mi = kgH2O ? amount/kgH2O : 0.0;
            return units::convert(mi, "molal", info.units);
        }
        if(info.quantity == "molarity")
        {
        	const Index iAqueous = system.indexPhaseWithError("Aqueous");
        	double amount = 0.0;
        	if(info.species.size())
				amount = state.speciesAmount(info.species);
        	if(info.element.size())
				amount = state.elementAmountInPhase(info.element, "Aqueous");
            const double volume = properties.phaseVolumes()[iAqueous].val;
            const double liter = convertCubicMeterToLiter(volume);
            const double ci = liter ? amount/liter : 0.0;
            return units::convert(ci, "molar", info.units);
        }
        if(info.quantity == "activity")
        {
            Index index = system.indexSpeciesWithError(info.species);
            const double ln_ai = properties.lnActivities().val[index];
            return std::exp(ln_ai);
        }
        if(info.quantity == "activityCoefficient")
        {
            Index index = system.indexSpeciesWithError(info.species);
            const double ln_gi = properties.lnActivityCoefficients().val[index];
            return std::exp(ln_gi);
        }
        if(info.quantity == "pH")
        {
            const Index iH = system.indexSpeciesWithError("H+");
            const double ln_aH = properties.lnActivities().val[iH];
            return -ln_aH/std::log(10);
        }
        if(info.quantity == "rate")
        {
            // Return zero if there are no reactions
            if(r.val.rows() == 0) return 0.0;
            Index index = reactions.indexReactionWithError(info.reaction);
            const double ri = r.val[index];
            return units::convert(ri, "mol/s", info.units);
        }

        RuntimeError("Cannot calculate the chemical quantity `" + str + "`.",
            "The string `" + str + "` does not contain a valid quantity.");

        return 0.0;
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

auto ChemicalQuantity::operator[](std::string quantity) const -> double
{
    return value(quantity);
}

} // namespace Reaktoro
