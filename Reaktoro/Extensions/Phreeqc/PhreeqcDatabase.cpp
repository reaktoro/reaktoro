// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "PhreeqcDatabase.hpp"

// C++ includes
#include <iostream>
#include <iomanip>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcUtils.hpp>

namespace Reaktoro {

struct PhreeqcDatabase::Impl
{
	/// The PHREEQC instance
	PHREEQC phreeqc;

	/// The list of elements in the database
	ElementList elements;

	/// The list of species in the database
	SpeciesList species;

	/// Construct a default PhreeqcDatabase::Impl object.
	Impl()
	{}

	/// Construct a PhreeqcDatabase::Impl object with given database.
	auto load(String database) -> void
	{
		// Load the PHREEQC database
		PhreeqcUtils::load(phreeqc, database);

		// Create the Element objects
		for(auto i = 0; i < phreeqc.count_elements; ++i)
			elements.append(createElement(phreeqc.elements[i]));

		// Create the Species objects using PhreeqcSpecies pointers (aqueous, exchange, surface species)
		for(auto i = 0; i < phreeqc.count_s; ++i)
			species.append(createSpecies(phreeqc.s[i]));

		// Create the Species objects using PhreeqcPhase pointers (gases and minerals)
		for(auto i = 0; i < phreeqc.count_phases; ++i)
			species.append(createSpecies(phreeqc.phases[i]));

		auto typestr = [](int val)
		{
			if(val == AQ       ) return "AQ";
			if(val == HPLUS    ) return "HPLUS";
			if(val == H2O      ) return "H2O";
			if(val == EMINUS   ) return "EMINUS";
			if(val == SOLID    ) return "SOLID";
			if(val == EX       ) return "EX";
			if(val == SURF     ) return "SURF";
			if(val == SURF_PSI ) return "SURF_PSI";
			if(val == SURF_PSI1) return "SURF_PSI1";
			if(val == SURF_PSI2) return "SURF_PSI2";
		};

		std::cout << "SPECIES" << '\n';
		for(auto i = 0; i < phreeqc.count_s; ++i)
		{
			std::cout << std::left << std::setw(20) << typestr(phreeqc.s[i]->type);
			std::cout << std::left << std::setw(20) << phreeqc.s[i]->name;
			std::cout << std::left << std::setw(20) << phreeqc.s[i]->z;
			std::cout << std::left << std::setw(20) << phreeqc.s[i]->logk[logK_T0];
			std::cout << std::left << std::setw(20) << phreeqc.s[i]->logk[delta_h];
			std::cout << std::left << std::setw(20) << phreeqc.s[i]->dha;
			std::cout << std::left << std::setw(20) << phreeqc.s[i]->dhb;
			std::cout << std::left << std::setw(20) << phreeqc.s[i]->a_f;
			std::cout << '\n';
		}

		std::cout << "PHASES" << '\n';
		for(auto i = 0; i < phreeqc.count_phases; ++i)
		{
			std::cout << std::left << std::setw(20) << typestr(phreeqc.phases[i]->type);
			std::cout << std::left << std::setw(20) << phreeqc.phases[i]->name;
			std::cout << std::left << std::setw(20) << phreeqc.phases[i]->formula;
			std::cout << std::left << std::setw(20) << phreeqc.phases[i]->logk[logK_T0];
			std::cout << std::left << std::setw(20) << phreeqc.phases[i]->logk[delta_h];
			std::cout << std::left << std::setw(20) << phreeqc.phases[i]->t_c;
			std::cout << '\n';
		}
	}

	/// Create an Element object with given PhreeqcElement pointer.
	auto createElement(const PhreeqcElement* e) const -> Element
	{
		Element element;
		element = element.withName(e->master->s->name);
		element = element.withSymbol(e->name);
		element = element.withMolarMass(e->gfw / 1000.0); // convert from g/mol to kg/mol
		return element;
	}

	/// Create the elemental composition of a Species object with given Phreeqc species pointer
	template<typename SpeciesType>
	auto createElementalComposition(const SpeciesType* s) const -> Map<Element, double>
	{
		Map<Element, double> pairs;
		for(auto&& [symbol, coeff] : PhreeqcUtils::elements(s))
			pairs.emplace(elements.get(symbol), coeff);
		return pairs;
	}

	/// Create the reactant species and their stoichiometries in a formation reaction.
	template<typename SpeciesType>
	auto createReactants(const SpeciesType* s) const -> Pairs<Species, double>
	{
		Pairs<Species, double> pairs;
		for(const auto& [reactant, coeff] : PhreeqcUtils::reactants(s))
			pairs.emplace_back(species.get(reactant), coeff);
		return pairs;
	}

	/// Create the equilibrium constant function (log base 10) of the formation reaction.
	template<typename SpeciesType>
	auto createEquilibriumConstantFn(const SpeciesType* s) const -> Fn<real(real, real)>
	{
		const auto fn = [=](real T, real P)
		{
			return PhreeqcUtils::lgEquilibriumConstant(s, T, P);
		};
		return fn;
	}

	/// Create the enthalpy change function of the formation reaction (in J/mol).
	template<typename SpeciesType>
	auto createEnthalpyChangeFn(const SpeciesType* s) const -> Fn<real(real, real)>
	{
		const auto fn = [=](real T, real P)
		{
			return PhreeqcUtils::enthalpyChange(s, T, P);
		};
		return fn;
	}

	/// Create the formation reaction of a product species.
	template<typename SpeciesType>
	auto createFormationReaction(const SpeciesType* s) const -> FormationReaction
	{
		FormationReaction reaction;
		reaction = reaction.withProduct(PhreeqcUtils::name(s));
		reaction = reaction.withReactants(createReactants(s));
		reaction = reaction.withEquilibriumConstantFn(createEquilibriumConstantFn(s));
		reaction = reaction.withEnthalpyChangeFn(createEnthalpyChangeFn(s));
		return reaction;
	}

	/// Create a Species object with given Phreeqc species pointer.
	template<typename SpeciesType>
	auto createSpecies(const SpeciesType* s) const -> Species
	{
		Species species;
		species = species.withName(PhreeqcUtils::name(s));
		species = species.withFormula(PhreeqcUtils::formula(s));
		species = species.withElements(createElementalComposition(s));
		species = species.withCharge(PhreeqcUtils::charge(s));
		species = species.withAggregateState(PhreeqcUtils::aggregateState(s));
		species = species.withFormationReaction(createFormationReaction(s));
		species = species.withAttachedData(s);
		return species;
	}
};

PhreeqcDatabase::PhreeqcDatabase()
: pimpl(new Impl())
{}

PhreeqcDatabase::PhreeqcDatabase(String database)
: PhreeqcDatabase()
{
	pimpl->load(database);
}

auto PhreeqcDatabase::elements() const -> ElementListConstRef
{
	return pimpl->elements;
}

auto PhreeqcDatabase::species() const -> SpeciesListConstRef
{
	return pimpl->species;
}

} // namespace Reaktoro
