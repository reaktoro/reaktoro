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

#include "PhreeqcUtils.hpp"

#ifdef LINK_PHREEQC

// Phreeqc includes
#define Phreeqc PHREEQC
#define protected public
#include <phreeqc/Phreeqc.h>
#include <phreeqc/GasPhase.h>
#undef Phreeqc

// Reaktoro includes
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Thermodynamics/Core/Database.hpp>
#include <Reaktoro/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/GaseousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/MineralSpecies.hpp>

namespace Reaktoro {

const double gram_to_kilogram = 1e-3;
const double kilogram_to_gram = 1e+3;
const double pascal_to_atm = 9.86923267e-6;
const double atm_to_pascal = 1/pascal_to_atm;

auto loadDatabase(PHREEQC& phreeqc, std::string filename) -> void
{
    std::istream* db_cookie = new std::ifstream(filename, std::ios_base::in);
    int errors = phreeqc.do_initialize();
    if(errors != 0) throw std::runtime_error("Error after call to `do_initialize` from Phreeqc instance.");
    phreeqc.Get_phrq_io()->push_istream(db_cookie);
    errors = phreeqc.read_database();
    phreeqc.Get_phrq_io()->clear_istream();
    phreeqc.Get_phrq_io()->Set_error_ostream(&std::cerr);
    phreeqc.Get_phrq_io()->Set_output_ostream(&std::cout);
    if(errors != 0) throw std::runtime_error("Error loading the database file in the Phreeqc instance.");
}

auto loadScript(PHREEQC& phreeqc, std::string filename) -> void
{
    std::istream* in_cookie = new std::ifstream(filename, std::ios_base::in);
    phreeqc.Get_phrq_io()->push_istream(in_cookie);
    int errors = phreeqc.run_simulations();
    phreeqc.Get_phrq_io()->clear_istream();
    if(errors != 0) throw std::runtime_error("Error when running the simulation with the Phreeqc instance.");
}

auto findElement(PHREEQC& phreeqc, std::string name) -> element*
{
    for(int i = 0; i < phreeqc.count_elements; ++i)
        if(phreeqc.elements[i]->name == name)
            return phreeqc.elements[i];
    return nullptr;
}

auto findSpecies(PHREEQC& phreeqc, std::string name) -> species*
{
    return phreeqc.s_search(name.c_str());
}

auto findPhase(PHREEQC& phreeqc, std::string name) -> phase*
{
    for(int i = 0; i < phreeqc.count_phases; ++i)
        if(phreeqc.phases[i]->name == name)
            return phreeqc.phases[i];
    return nullptr;
}

auto isGaseousSpecies(const phase* p) -> bool
{
	std::string name(p->name);
	return name.find("(g)") < name.size();
}

auto isMineralSpecies(const phase* p) -> bool
{
	return !isGaseousSpecies(p);
}

auto getElements(const species* s) -> std::map<element*, double>
{
    std::map<element*, double> elements;
    for(auto iter = s->next_elt; iter->elt != NULL; iter++)
        elements.emplace(iter->elt, iter->coef);
    return elements;
}

auto getElements(const phase* p) -> std::map<element*, double>
{
    std::map<element*, double> elements;
    for(auto iter = p->next_elt; iter->elt != NULL; iter++)
        elements.emplace(iter->elt, iter->coef);
    return elements;
}

auto getElementStoichiometry(std::string element, const species* s) -> double
{
	// Return species charge if element is electrical charge Z
    if(element == "Z") return s->z;

    // Iterate over all elements in the species
    for(auto iter = s->next_elt; iter->elt != nullptr; iter++)
        if(iter->elt->name == element)
            return iter->coef;

    return 0.0;
}

auto getElementStoichiometry(std::string element, const phase* p) -> double
{
    // Iterate over all elements in the phase (gasesous or mineral species)
    for(auto iter = p->next_elt; iter->elt != nullptr; iter++)
        if(iter->elt->name == element)
            return iter->coef;
    return 0.0;
}

auto getReactionEquation(const species* s) -> std::map<std::string, double>
{
    std::map<std::string, double> equation;

    // Check if there is any reaction defined by this species
    if(s->rxn_x)
    {
        // Iterate over all reactants
        for(auto iter = s->rxn_x->token; iter->name != nullptr; iter++)
            equation.emplace(iter->name, iter->coef);

        // Check if the reaction is defined by only one species
        if(equation.size() <= 1)
            equation.clear();
    }

    return equation;
}

auto getReactionEquation(const phase* p) -> std::map<std::string, double>
{
    std::map<std::string, double> equation;
    for(auto iter = p->rxn_x->token; iter->name != nullptr; iter++)
        equation.emplace(iter->name, iter->coef);
    return equation;
}

auto collectAqueousSpecies(PHREEQC& phreeqc) -> std::vector<species*>
{
    std::set<species*> aqueous_species;
    for(int i = 0; i < phreeqc.count_s; ++i)
    {
        // Check if the species is considered in the model
        if(phreeqc.s[i]->in)
        {
            aqueous_species.insert(phreeqc.s[i]);

            // Add the reactants in the reaction defined by this species
            for(auto pair : getReactionEquation(phreeqc.s[i]))
                aqueous_species.insert(findSpecies(phreeqc, pair.first));
        }
    }

    return {aqueous_species.begin(), aqueous_species.end()};
}

auto collectSecondarySpecies(PHREEQC& phreeqc) -> std::vector<species*>
{
    std::vector<species*> secondary_species;
    for(species* s : collectAqueousSpecies(phreeqc))
        if(getReactionEquation(s).size())
            secondary_species.push_back(s);

    return secondary_species;
}

auto collectGaseousSpecies(PHREEQC& phreeqc) -> std::vector<phase*>
{
    std::set<phase*> gaseous_species;

    // Collect the gaseous species of the cxxGasPhase instances
    for(const auto& pair : phreeqc.Rxn_gas_phase_map)
    {
        const cxxGasPhase& gas_phase = pair.second;
        for(const cxxGasComp& component : gas_phase.Get_gas_comps())
            gaseous_species.insert(findPhase(phreeqc, component.Get_phase_name()));
    }

    return {gaseous_species.begin(), gaseous_species.end()};
}

auto collectMineralSpecies(PHREEQC& phreeqc) -> std::vector<phase*>
{
    std::vector<phase*> pure_minerals;

    unknown** x = phreeqc.x;
    for(int i = 0; i < phreeqc.count_unknowns; ++i)
    	if(x[i]->type == PP && x[i]->phase->rxn_x != nullptr && x[i]->phase->in == true)
    		pure_minerals.push_back(x[i]->phase);

    return pure_minerals;
}

auto collectGaseousSpeciesInSpeciationList(PHREEQC& phreeqc) -> std::vector<phase*>
{
    std::vector<phase*> gaseous_species;

    auto is_gas = [](phase* p)
    {
        std::string name(p->name);
        const bool is_in = p->in;
        const bool is_gas = name.find("(g)") < name.size();
        return is_in && is_gas;
    };

    for(int i = 0; i < phreeqc.count_phases; ++i)
        if(is_gas(phreeqc.phases[i]))
            gaseous_species.push_back(phreeqc.phases[i]);

    return gaseous_species;
}

auto collectMineralSpeciesInSpeciationList(PHREEQC& phreeqc) -> std::vector<phase*>
{
    auto is_mineral = [](phase* p)
    {
        std::string name(p->name);
        const bool is_in = p->in;
        const bool is_solid = p->type == SOLID;
        const bool is_not_gas = name.find("(g)") > name.size();
        return is_in && is_solid && is_not_gas;
    };

    std::vector<phase*> pure_minerals;
    for(int i = 0; i < phreeqc.count_phases; ++i)
        if(is_mineral(phreeqc.phases[i]))
            pure_minerals.push_back(phreeqc.phases[i]);
    return pure_minerals;
}

auto index(std::string name, const std::vector<species*>& pointers) -> unsigned
{
    auto compare = [=](species* entry) { return entry->name == name; };
    return std::find_if(pointers.begin(), pointers.end(), compare) - pointers.begin();
}

auto index(std::string name, const std::vector<phase*>& pointers) -> unsigned
{
    auto compare = [=](phase* entry) { return entry->name == name; };
    return std::find_if(pointers.begin(), pointers.end(), compare) - pointers.begin();
}

auto speciesAmountsInSpecies(const std::vector<species*>& pointers) -> Vector
{
    const unsigned size = pointers.size();
    Vector n(size);
    for(unsigned i = 0; i < size; ++i)
        n[i] = pointers[i]->moles;
    return n;
}

auto speciesAmountsInPhases(const std::vector<phase*>& pointers) -> Vector
{
    const unsigned size = pointers.size();
    Vector n(size);
    for(unsigned i = 0; i < size; ++i)
        n[i] = pointers[i]->moles_x;
    return n;
}

auto lnEquilibriumConstant(const double* l_logk, double T, double P) -> double
{
    // The universal gas constant defined in Phreeqc (in units of kJ/(K*mol))
    const double R = R_KJ_DEG_MOL;
    const double ln_10 = std::log(10.0);

    // Calculate log10(k) for this temperature and pressure
    double log10_k = l_logk[logK_T0]
        - l_logk[delta_h] * (298.15 - T)/(ln_10*R*T*298.15)
        + l_logk[T_A1]
        + l_logk[T_A2]*T
        + l_logk[T_A3]/T
        + l_logk[T_A4]*std::log10(T)
        + l_logk[T_A5]/(T*T)
        + l_logk[T_A6]*(T*T);

    return log10_k * ln_10;
}

auto lnEquilibriumConstant(const species* s, double T, double P) -> double
{
    return lnEquilibriumConstant(s->rxn_x->logk, T, P);
}

auto lnEquilibriumConstant(const phase* p, double T, double P) -> double
{
    return lnEquilibriumConstant(p->rxn_x->logk, T, P);
}

auto convert(const PHREEQC& phreeqc, const Database& master_database, const std::vector<double>& tpoints, const std::vector<double>& ppoints) -> Database
{
	auto reaktoro_species_name = [](const species* s) -> std::string
	{
		const std::string name = s->name;
		const double charge = s->z;
		if(name == "H2O") return "H2O(l)";
		if(charge == 0) return name + "(aq)";
		if(charge < 0) return name.substr(0, name.rfind('-')) + std::string('-', std::abs(charge));
		else return name.substr(0, name.rfind('+')) + std::string('+', std::abs(charge));
	};

	auto index_master_species = [&](std::string name)
	{
		for(int i = 0; i < phreeqc.count_master; ++i)
			if(phreeqc.master[i]->s->name == name)
				return i;
		return phreeqc.count_master;
	};

	auto is_master_species = [&](const species* s)
	{
		return index_master_species(s->name) < phreeqc.count_master;
	};

	auto create_element = [&](const element* e) -> Element
	{
		Element element;
		element.setName(e->name);
		element.setMolarMass(e->gfw);
		return element;
	};

	auto elements_in_species = [&](const species* s)
	{
		std::map<Element, double> res;
		for(auto pair : getElements(s))
			res.insert({create_element(pair.first), pair.second});
		return res;
	};

	auto elements_in_phase = [&](const phase* p)
	{
		std::map<Element, double> res;
		for(auto pair : getElements(p))
			res.insert({create_element(pair.first), pair.second});
		return res;
	};

	auto thermo_data_in_product_aqueous_species = [&](const species* s)
	{
		auto lnk_fn = [&](double T, double P)
		{
			return -lnEquilibriumConstant(s, T, P);
		};

		ReactionThermoInterpolatedProperties rprops;
		rprops.equation = ReactionEquation(getReactionEquation(s));
		rprops.lnk = BilinearInterpolator(tpoints, ppoints, lnk_fn);

		AqueousSpeciesThermoData data;
		data.reaction.set(rprops);

		return data;
	};

	auto thermo_data_in_gaseous_species = [&](const phase* p)
	{
		auto lnk_fn = [&](double T, double P)
		{
			return lnEquilibriumConstant(p, T, P);
		};

		ReactionThermoInterpolatedProperties rprops;
		rprops.equation = ReactionEquation(getReactionEquation(p));
		rprops.lnk = BilinearInterpolator(tpoints, ppoints, lnk_fn);

		GaseousSpeciesThermoData data;
		data.reaction.set(rprops);

		return data;
	};

	auto thermo_data_in_mineral_species = [&](const phase* p)
	{
		auto lnk_fn = [&](double T, double P)
		{
			return lnEquilibriumConstant(p, T, P);
		};

		const double cm3_to_m3 = 1e-6;
		const double molar_volume = p->logk[vm0] * cm3_to_m3;

		SpeciesThermoInterpolatedProperties sprops;
		sprops.volume = BilinearInterpolator({278.15}, {1e5}, {molar_volume});

		ReactionThermoInterpolatedProperties rprops;
		rprops.equation = ReactionEquation(getReactionEquation(p));
		rprops.lnk = BilinearInterpolator(tpoints, ppoints, lnk_fn);

		MineralSpeciesThermoData data;
		data.properties.set(sprops);
		data.reaction.set(rprops);

		return data;
	};

	auto create_aqueous_species = [&](const species* s) -> AqueousSpecies
	{
		const std::string name = reaktoro_species_name(s);

		if(is_master_species(s))
			return master_database.aqueousSpecies(name);
		else
		{
			AqueousSpecies species;
			species.setName(name);
			species.setCharge(s->z);
			species.setElements(elements_in_species(s));
			species.setThermoData(thermo_data_in_product_aqueous_species(s));
			return species;
		}
	};

	auto create_gaseous_species = [&](const phase* p) -> GaseousSpecies
	{
		GaseousSpecies species;
		species.setName(p->name);
		species.setElements(elements_in_phase(p));
		species.setThermoData(thermo_data_in_gaseous_species(p));
		return species;
	};

	auto create_mineral_species = [&](const phase* p) -> MineralSpecies
	{
		MineralSpecies species;
		species.setName(p->name);
		species.setElements(elements_in_phase(p));
		species.setThermoData(thermo_data_in_mineral_species(p));
		return species;
	};

	// Create the Database instance
	auto create_database = [&]()
	{
		Database database;

		for(int i = 0; i < phreeqc.count_s; ++i)
			database.addAqueousSpecies(create_aqueous_species(phreeqc.s[i]));

		for(int i = 0; i < phreeqc.count_phases; ++i)
			if(isGaseousSpecies(phreeqc.phases[i]))
				database.addGaseousSpecies(create_gaseous_species(phreeqc.phases[i]));
			else
				database.addMineralSpecies(create_mineral_species(phreeqc.phases[i]));

		return database;
	};

	return create_database();
}

} // namespace Reaktoro

#endif
