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

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/AggregateState.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Species.hpp>

// Phreeqc includes
#include <Reaktoro/Interfaces/PhreeqcLegacy.hpp>
#include <Reaktoro/Interfaces/PhreeqcUtils.hpp>

namespace Reaktoro {
namespace {

/// Convert a PhreeqcElement pointer into a Reaktoro::Element object.
auto convert(const PhreeqcElement* e) -> Element
{
	Element element;
	element = element.withName(e->name);
	element = element.withSymbol(e->name);
	element = element.withMolarMass(e->gfw);
	return element;
}

/// Create the Element objects and their coefficients with given PhreeqcSpecies pointer.
auto createElements(const PhreeqcSpecies* s) -> Species::Elements
{
	Species::Elements elements;
	for(auto&& [element, coeff] : PhreeqcUtils::elements(s))
		elements.emplace(convert(element), coeff);
	return elements;
}

/// Create the Element objects and their coefficients with given PhreeqcPhase pointer.
auto createElements(const PhreeqcPhase* p) -> Species::Elements
{
	Species::Elements elements;
	for(auto&& [element, coeff] : PhreeqcUtils::elements(p))
		elements.emplace(convert(element), coeff);
	return elements;
}

template<typename SpeciesType>
auto speciesThermoParamsPhreeqc(const SpeciesType& species) -> SpeciesThermoParamsPhreeqc
{
	SpeciesThermoParamsPhreeqc params;
	params.reaction.equation = PhreeqcUtils::reactionEquation(species);
	params.reaction.log_k = species->logk[logK_T0];
	params.reaction.delta_h = species->logk[delta_h];
	params.reaction.analytic = {
		species->logk[T_A1],
		species->logk[T_A2],
		species->logk[T_A3],
		species->logk[T_A4],
		species->logk[T_A5],
		species->logk[T_A6]};
	return params;
}

auto createAqueousSpecies(const PhreeqcSpecies* s) -> Species
{
	Species species;
	species = species.withName(s->name);
	species = species.withFormula(s->name);
	species = species.withElements(createElements(s));
	species = species.withCharge(s->z);
	species = species.withAggregateState(AggregateState::Aqueous);
	species = species.withAttachedData(speciesThermoParamsPhreeqc(s));
	return species;
}

auto createGaseousSpecies(const PhreeqcPhase* p) -> Species
{
	Species species;
	species = species.withName(p->name);
	species = species.withFormula(p->formula);
	species = species.withElements(createElements(p));
	species = species.withCharge(0.0);
	species = species.withAggregateState(AggregateState::Gas);
	species = species.withAttachedData(speciesThermoParamsPhreeqc(p));
	return species;
}

auto createMineralSpecies(const PhreeqcPhase* p) -> Species
{
	Species species;
	species = species.withName(p->name);
	species = species.withFormula(p->formula);
	species = species.withElements(createElements(p));
	species = species.withCharge(0.0);
	species = species.withAggregateState(AggregateState::Solid);
	species = species.withTags({"mineral"});
	species = species.withAttachedData(speciesThermoParamsPhreeqc(p));
	return species;
}

} // namespace

struct PhreeqcDatabase::Impl
{
	/// The PHREEQC instance
	PHREEQC phreeqc;

	/// The list of elements in the database
	std::vector<Element> elements;

	/// The list of aqueous species in the database
	std::vector<Species> aqueous_species;

	/// The list of gaseous species in the database
	std::vector<Species> gaseous_species;

	/// The list of mineral species in the database
	std::vector<Species> mineral_species;

	/// The aqueous master species in the database
	std::set<std::string> master_species;

	/// The indices of the aqueous master species in the database
	Indices idx_master_species;

	/// The map from a master species to the product species composed by it
	std::vector<std::set<std::string>> from_master_to_product_species;

	auto load(std::string filename) -> void
	{
		// Clear current state
		elements.clear();
		aqueous_species.clear();
		gaseous_species.clear();
		mineral_species.clear();
		master_species.clear();
		idx_master_species.clear();
		from_master_to_product_species.clear();

		// Initialize the phreeqc instance
		int errors = phreeqc.do_initialize();
		Assert(errors == 0, "Could not initialize Phreeqc.",
			"Call to method `Phreeqc::do_initialize` failed.");

		// Initialize the phreeqc database
		std::istream* db_cookie = new std::ifstream(filename, std::ios_base::in);
		phreeqc.Get_phrq_io()->push_istream(db_cookie);
		errors = phreeqc.read_database();
		phreeqc.Get_phrq_io()->clear_istream();
		phreeqc.Get_phrq_io()->Set_error_ostream(&std::cerr);
		phreeqc.Get_phrq_io()->Set_output_ostream(&std::cout);
		Assert(errors == 0, "Could not load the Phreeqc database file `" + filename + "`.",
			"Ensure `" + filename + "` points to the right path to the database file.");

		// Initialize the set of master species
		for(int i = 0; i < phreeqc.count_master; ++i)
			master_species.insert(phreeqc.master[i]->s->name);

		// Initialize the indices of the aqueous species that are master species
		for(int i = 0; i < phreeqc.count_s; ++i)
			if(master_species.count(phreeqc.s[i]->name))
				idx_master_species.push_back(i);

		// Initialize the map from master species to the product species
		from_master_to_product_species.resize(master_species.size());

		for(unsigned i = 0; i < idx_master_species.size(); ++i)
		{
		    const Index ispecies = idx_master_species[i];
			const std::string master_name = phreeqc.s[ispecies]->name;

			// Loop over all Phreeqc species instances (aqueous species)
			for(int j = 0; j < phreeqc.count_s; ++j)
			{
				const std::string product_name = phreeqc.s[j]->name;
				const auto equation = PhreeqcUtils::reactionEquation(phreeqc.s[j]);
				for(auto pair : equation)
					if(pair.first == master_name)
						from_master_to_product_species[i].insert(product_name);
			}

			// Loop over all Phreeqc phase instances (gaseous or mineral species)
			for(int j = 0; j < phreeqc.count_phases; ++j)
			{
				const std::string product_name = phreeqc.phases[j]->name;
				const auto equation = PhreeqcUtils::reactionEquation(phreeqc.phases[j]);
				for(auto pair : equation)
					if(pair.first == master_name)
						from_master_to_product_species[i].insert(product_name);
			}
		}

		// Initialize the elements
		for(int i = 0; i < phreeqc.count_elements; ++i)
			elements.push_back(convert(phreeqc.elements[i]));

		// Initialize the aqueous species
		for(int i = 0; i < phreeqc.count_s; ++i)
			aqueous_species.push_back(createAqueousSpecies(phreeqc.s[i]));

		// Initialize the gaseous and mineral species
		for(int i = 0; i < phreeqc.count_phases; ++i)
			if(PhreeqcUtils::isGaseousSpecies(phreeqc.phases[i]))
				gaseous_species.push_back(createGaseousSpecies(phreeqc.phases[i]));
			else
				mineral_species.push_back(createMineralSpecies(phreeqc.phases[i]));
	}
};

PhreeqcDatabase::PhreeqcDatabase()
: pimpl(new Impl())
{}

PhreeqcDatabase::PhreeqcDatabase(std::string filename)
: PhreeqcDatabase()
{
	load(filename);
}

auto PhreeqcDatabase::load(std::string filename) -> void
{
	pimpl->load(filename);
}

auto PhreeqcDatabase::masterSpecies() const -> std::set<std::string>
{
	return pimpl->master_species;
}

auto PhreeqcDatabase::aqueousSpecies() const -> const std::vector<Species>&
{
	return pimpl->aqueous_species;
}

auto PhreeqcDatabase::gaseousSpecies() const -> const std::vector<Species>&
{
	return pimpl->gaseous_species;
}

auto PhreeqcDatabase::mineralSpecies() const -> const std::vector<Species>&
{
	return pimpl->mineral_species;
}

} // namespace Reaktoro
