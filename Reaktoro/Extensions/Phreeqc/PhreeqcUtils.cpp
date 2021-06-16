// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

#include "PhreeqcUtils.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Constants.hpp>

namespace Reaktoro {
namespace PhreeqcUtils {

auto load(PHREEQC& phreeqc, String database) -> void
{
    // Initialize the phreeqc instance
    int errors = phreeqc.do_initialize();
    Assert(errors == 0, "Could not initialize PHREEQC.",
        "Call to method `Phreeqc::do_initialize` failed.");

    //-------------------------------------------------------------------------
    // Check below if the argument `database` is:
    // 1) a path to a database file; or
    // 2) a multi-line string containing the database contents it self.
    //-------------------------------------------------------------------------

    // Initialize the phreeqc database
    std::istream* db_cookie;
    if(database.find('\n') == String::npos) // check if there no newline char in argument `database`
        db_cookie = new std::ifstream(database, std::ios_base::in); // the database argument is a path to a file
    else db_cookie = new std::istringstream(database); // the database argument contains the database contents itself

    phreeqc.Get_phrq_io()->push_istream(db_cookie);
    errors = phreeqc.read_database();
    phreeqc.Get_phrq_io()->clear_istream();
    phreeqc.Get_phrq_io()->Set_error_ostream(&std::cerr);
    phreeqc.Get_phrq_io()->Set_output_ostream(&std::cout);
    Assert(errors == 0, "Could not load the PHREEQC database file `" + database + "`.",
        "Ensure `" + database + "` points to the right path to the database file.");
}

auto execute(PHREEQC& phreeqc, String input, String output) -> void
{
    //-------------------------------------------------------------------------
    // Check below if the argument `input` is:
    // 1) a path to a input script file; or
    // 2) a multi-line string containing the input script contents it self.
    //-------------------------------------------------------------------------
    std::istream* in_cookie;
    if(input.find('\n') == String::npos) // check if there no newline char in argument `input`
        in_cookie = new std::ifstream(input, std::ios_base::in); // the input argument is a path to a file
    else in_cookie = new std::istringstream(input); // the input argument contains the input script contents itself

    // Check if output to a file should be performed
    std::ofstream* out = output.empty() ? nullptr : new std::ofstream(output);

    // Set the output and error streams
    phreeqc.Get_phrq_io()->Set_output_ostream(out);
    phreeqc.Get_phrq_io()->Set_error_ostream(out);

    // Set the input stream and execute the simulation
    phreeqc.Get_phrq_io()->push_istream(in_cookie);
    int errors = phreeqc.run_simulations();
    phreeqc.Get_phrq_io()->clear_istream();

    // Delete the dynamically allocated stream objects
    if(out != nullptr) delete out;

    // Check if there was any error
    Assert(errors == 0, "Failed to execute the PHREEQC input script `" + input + "`.",
        "There was a Phreeqc error when executing this input script file.");
}

auto findElement(const PHREEQC& phreeqc, String name) -> PhreeqcElement*
{
    for(auto i = 0; i < phreeqc.count_elements; ++i)
        if(phreeqc.elements[i]->name == name)
            return phreeqc.elements[i];
    return nullptr;
}

auto findSpecies(const PHREEQC& phreeqc, String name) -> PhreeqcSpecies*
{
    for(auto i = 0; i < phreeqc.count_s; ++i)
        if(phreeqc.s[i]->name == name)
            return phreeqc.s[i];
    return nullptr;
}

auto findPhase(const PHREEQC& phreeqc, String name) -> PhreeqcPhase*
{
    for(auto i = 0; i < phreeqc.count_phases; ++i)
        if(phreeqc.phases[i]->name == name)
            return phreeqc.phases[i];
    return nullptr;
}

auto isAqueousSpecies(const PhreeqcSpecies* species) -> bool
{
    return oneof(species->type, AQ, HPLUS, H2O, EMINUS);
}

auto isGaseousSpecies(const PhreeqcPhase* phase) -> bool
{
    String name(phase->name);
    return name.find("(g)") < name.size();
}

auto isMineralSpecies(const PhreeqcPhase* phase) -> bool
{
    return !isGaseousSpecies(phase) && phase->type == SOLID;

    // Note 1: Gases are marked SOLID as well. That's why we also check whether
    // this phase pointer is not a gaseous species!

    // Note 2: Exchange species are also present in `PHREEQC::phases` container.
    // These are marked as EX. Thus, checking only if not a gaseous species is
    // not sufficient to determine if the phase is a mineral species.
}

auto isExchangeSpecies(const PhreeqcSpecies* species) -> bool
{
    return species->type == EX;
}

auto isSurfaceSpecies(const PhreeqcSpecies* species) -> bool
{
    return oneof(species->type, SURF, SURF_PSI, SURF_PSI1, SURF_PSI2);
}

auto symbol(const PhreeqcElement* element) -> String
{
    return element->name; // element symbol is H, H(0), C, Fe, Fe(3), etc.
}

auto name(const PhreeqcElement* element) -> String
{
    return element->master->s->name; // element name is H+, H2, CO3-2, Fe+2, Fe+3, etc.
}

auto molarMass(const PhreeqcElement* element) -> double
{
    return element->gfw / 1000.0; // convert molar mass from g/mol to kg/mol
}

auto name(const PhreeqcSpecies* species) -> String
{
    return species->name;
}

auto name(const PhreeqcPhase* phase) -> String
{
    return phase->name;
}

auto formula(const PhreeqcSpecies* species) -> String
{
    return species->name;
}

auto formula(const PhreeqcPhase* phase) -> String
{
    return phase->formula;
}

auto elements(const PhreeqcSpecies* species) -> Map<PhreeqcElement*, double>
{
    Map<PhreeqcElement*, double> elements;
    for(auto iter = species->next_elt; iter->elt != nullptr; iter++)
        elements.emplace(iter->elt, iter->coef);
    return elements;
}

auto elements(const PhreeqcPhase* phase) -> Map<PhreeqcElement*, double>
{
    Map<PhreeqcElement*, double> elements;
    for(auto iter = phase->next_elt; iter->elt != nullptr; iter++)
        elements.emplace(iter->elt, iter->coef);
    return elements;
}

auto charge(const PhreeqcSpecies* species) -> double
{
    return species->z;
}

auto charge(const PhreeqcPhase* phase) -> double
{
    return 0.0;
}

auto aggregateState(const PhreeqcSpecies* species) -> AggregateState
{
    if(isAqueousSpecies(species))  return AggregateState::Aqueous;
    if(isExchangeSpecies(species)) return AggregateState::IonExchange;
    if(isSurfaceSpecies(species))  return AggregateState::Adsorbed;
    return AggregateState::Undefined;
}

auto aggregateState(const PhreeqcPhase* phase) -> AggregateState
{
    if(isGaseousSpecies(phase)) return AggregateState::Gas;
    if(isMineralSpecies(phase)) return AggregateState::Solid;
    return AggregateState::Undefined;
}

template<typename SpeciesType>
auto stoichiometryAux(const SpeciesType* species, String name) -> double
{
    const auto equation = reactionEquation(species);
    for(const auto& [x, y] : equation)
        if(x == name)
            return y;
    return 0.0;
}

auto stoichiometry(const PhreeqcSpecies* species, String name) -> double
{
    return stoichiometryAux(species, name);
}

auto stoichiometry(const PhreeqcPhase* phase, String name) -> double
{
    return stoichiometryAux(phase, name);
}

auto reactionEquation(const PhreeqcSpecies* species) -> Pairs<String, double>
{
    // Check if there is any reaction defined by this species.
    if(species->rxn_x == nullptr) return {};

    // The reaction equation as pairs of names and stoichiometries
    Pairs<String, double> pairs;

    // Iterate over all species in the reaction, get their names and stoichiometries.
    for(auto iter = species->rxn_x->token; iter->s != nullptr; iter++)
        pairs.emplace_back(iter->s->name, -iter->coef);

    // Note: The minus sign on `iter->coef` above is needed to reverse the
    // reaction so that method PhreeqcUtils::lgEquilibriumConstant returns a
    // consistent logK value.

    // Check if the reaction is a trivial reaction (e.g., H+ = H+)
    // Phreeqc generates this reactions for master species. Here, we
    // return a empty reaction equation in this case.
    if(pairs.size() <= 1) return {};

    return pairs;
}

auto reactionEquation(const PhreeqcPhase* phase) -> Pairs<String, double>
{
    // Check if there is any reaction defined by this species.
    if(phase->rxn_x == nullptr) return {};

    // The reaction equation as pairs of names and stoichiometries
    Pairs<String, double> pairs;

    // Iterate over all species in the reaction, get their names and stoichiometries.
    // Note that the reaction equation for a Phreeqc phase is read differently from a
    // Phreeqc species instance.
    pairs.emplace_back(phase->rxn_x->token->name, phase->rxn_x->token->coef);
    for(auto iter = phase->rxn_x->token + 1; iter->s != nullptr; iter++)
        pairs.emplace_back(iter->s->name, iter->coef);

    // Check if the reaction is a trivial reaction (e.g., X- = X-)
    // Phreeqc generates this reactions for master species. Here, we
    // return a empty reaction equation in this case.
    if(pairs.size() <= 1) return {};

    return pairs;
}

auto reactants(const PhreeqcSpecies* species) -> Pairs<PhreeqcSpecies*, double>
{
    // Assert there is a reaction defined for this species (even if a master species).
    assert(species->rxn != nullptr);

    // The reactants in the formation reaction of the species (and their stoichiometric coefficients)
    Pairs<PhreeqcSpecies*, double> reactants;

    // Iterate over all species in the reaction, get their names and stoichiometries.
    for(auto iter = species->rxn->token; iter->s != nullptr; iter++)
        if(iter->s->name != species->name)
            reactants.emplace_back(iter->s, iter->coef);

    // Note: The species in the reaction with same name as the product
    // species is skiped from the list of reactants.

    return reactants;
}

auto reactants(const PhreeqcPhase* phase) -> Pairs<PhreeqcSpecies*, double>
{
    // Assert there is a reaction defined for this species.
    assert(phase->rxn != nullptr);

    // The reactants in the formation reaction of the species (and their stoichiometric coefficients)
    Pairs<PhreeqcSpecies*, double> reactants;

    // Iterate over all species in the reaction, get their names and stoichiometries.
    for(auto iter = phase->rxn->token + 1; iter->s != nullptr; iter++)
        reactants.emplace_back(iter->s, iter->coef);

    // Note: The species in the reaction with same name as the product species
    // is skiped from the list of reactants.

    return reactants;
}

auto isMasterSpecies(const PhreeqcSpecies* species) -> bool
{
    return reactants(species).empty();
}

auto isMasterSpecies(const PhreeqcPhase* phase) -> bool
{
    return reactants(phase).empty();
}

auto index(String name, const Vec<PhreeqcSpecies*>& species) -> std::size_t
{
    auto compare = [=](PhreeqcSpecies* entry) { return entry->name == name; };
    return std::find_if(species.begin(), species.end(), compare) - species.begin();
}

auto index(String name, const Vec<PhreeqcPhase*>& phases) -> std::size_t
{
    auto compare = [=](PhreeqcPhase* entry) { return entry->name == name; };
    return std::find_if(phases.begin(), phases.end(), compare) - phases.begin();
}

auto activeAqueousSpecies(const PHREEQC& phreeqc) -> Vec<PhreeqcSpecies*>
{
    // Inspired by method `int Phreeqc::print_species(void)`
    std::set<PhreeqcSpecies*> species;

    // Loop over all species in `species_list` data-member of PHREEQC
    for(auto i = 0; i < phreeqc.count_species_list; ++i)
        if(phreeqc.species_list[i].s->type != EX &&
           phreeqc.species_list[i].s->type != SURF)
            species.insert(phreeqc.species_list[i].s);

    // Add electron species e- if present in some reaction
    for(auto s : species)
        if(stoichiometry(s, "e-") != 0.0)
            species.insert(phreeqc.s_eminus);

    return {species.begin(), species.end()};
}

auto activeExchangeSpecies(const PHREEQC& phreeqc) -> Vec<PhreeqcSpecies*>
{
    // Inspired by method `int Phreeqc::print_exchange(void)`
    std::set<PhreeqcSpecies*> species;
    for(auto i = 0; i < phreeqc.count_species_list; ++i)
        if(phreeqc.species_list[i].s->type == EX)
            species.insert(phreeqc.species_list[i].s);
    return {species.begin(), species.end()};
}

auto activeProductSpecies(const PHREEQC& phreeqc) -> Vec<PhreeqcSpecies*>
{
    Vec<PhreeqcSpecies*> species;

    // Loop over all species in `s_x` data-member of PHREEQC
    for(auto i = 0; i < phreeqc.count_s_x; ++i)
        if(phreeqc.s_x[i]->type != EX && phreeqc.s_x[i]->type != SURF) // exchange and surface species are currently not supported
            if(!reactionEquation(phreeqc.s_x[i]).empty()) //  a product species has a non-empty reaction equation
                species.push_back(phreeqc.s_x[i]);

    return species;
}

auto activeGaseousSpecies(const PHREEQC& phreeqc) -> Vec<PhreeqcPhase*>
{
    // Collect the gaseous species of the cxxGasPhase instances
    std::set<PhreeqcPhase*> gases;
    for(const auto& pair : phreeqc.Rxn_gas_phase_map)
    {
        const cxxGasPhase& gas_phase = pair.second;
        for(const cxxGasComp& component : gas_phase.Get_gas_comps())
            gases.insert(findPhase(phreeqc, component.Get_phase_name()));
    }
    return {gases.begin(), gases.end()};
}

auto activePhasesInEquilibriumPhases(const PHREEQC& phreeqc) -> Vec<PhreeqcPhase*>
{
    // Inspired by method `int Phreeqc::print_pp_assemblage(void)`
    Vec<PhreeqcPhase*> phases;
    unknown** x = phreeqc.x;
    for(auto i = 0; i < phreeqc.count_unknowns; ++i)
        if(x[i]->type == PP && x[i]->phase->rxn_x != nullptr && x[i]->phase->in != 0)
            phases.push_back(x[i]->phase);
    return phases;
}

auto activePhasesInSaturationList(const PHREEQC& phreeqc) -> Vec<PhreeqcPhase*>
{
    Vec<PhreeqcPhase*> phases;
    for(auto i = 0; i < phreeqc.count_phases; ++i)
        if(phreeqc.phases[i]->in)
            phases.push_back(phreeqc.phases[i]);
    return phases;
}

auto speciesAmounts(const PHREEQC& phreeqc, const Vec<PhreeqcSpecies*>& species) -> ArrayXd
{
    const auto num_species = species.size();
    ArrayXd n(num_species);
    for(auto i = 0; i < num_species; ++i)
        n[i] = species[i]->moles;
    return n;
}

auto speciesAmounts(const PHREEQC& phreeqc, const Vec<PhreeqcPhase*>& phases) -> ArrayXd
{
    const auto num_phases = phases.size();
    const auto num_unknowns = phreeqc.count_unknowns;
    ArrayXd n(num_phases);
    for(auto i = 0; i < num_phases; ++i)
    {
        n[i] = phases[i]->moles_x;

        // Check if the phase has zero molar amount
        if(n[i] <= 0.0)
        {
            // Check if there is any unknown for current phase and get its moles from there
            for(auto j = 0; j < num_unknowns; ++j)
            {
                if(phases[i] == phreeqc.x[j]->phase)
                {
                    n[i] = phreeqc.x[j]->moles;
                    break;
                }
            }
        }
    }

    return n;
}

} // namespace PhreeqcUtils
} // namespace Reaktoro
