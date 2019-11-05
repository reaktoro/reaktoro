// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include "Phreeqc.hpp"

// C++ includes
#include <map>

// Eigen includes
#include <Reaktoro/deps/eigen3/Eigen/Dense>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/ConvertUtils.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/ReactionEquation.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Interfaces/PhreeqcUtils.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

// Phreeqc includes
#define Phreeqc PHREEQC
#define protected public
#include <phreeqc/GasPhase.h>
#include <phreeqc/Phreeqc.h>
#undef Phreeqc

namespace Reaktoro {
namespace {

const double gram_to_kilogram = 1e-3;
const double pascal_to_bar = 1e-5;
const double pascal_to_atm = 9.86923267e-6;
const double atm_to_pascal = 1 / pascal_to_atm;
const double cm3_to_m3 = 1e-6;

// The critical properties of some gases
// The data are: {temperature [C], pressure [bar], acentric factor}
std::map<std::string, std::vector<double>> critical_properties =
    {
        {"Ar(g)", {150.9, 48.98, 0.0}},
        {"CH4(g)", {190.6, 45.99, 0.012}},
        {"C6H6O(g)", {694.3, 61.3, 0.444}},
        {"CO(g)", {132.9, 34.99, 0.048}},
        {"CO2(g)", {304.2, 73.83, 0.224}},
        {"C2H4(g)", {282.3, 50.4, 0.087}},
        {"H2(g)", {33.19, 13.13, -0.216}},
        {"H2O(g)", {647.1, 220.55, 0.345}},
        {"H2S(g)", {373.5, 89.63, 0.094}},
        {"He(g)", {5.2, 2.28, -0.39}},
        {"Kr(g)", {209.4, 55.02, 0.0}},
        {"N2(g)", {126.2, 34.0, 0.038}},
        {"Ne(g)", {44.0, 27.0, 0.0}},
        {"NH3(g)", {405.7, 112.8, 0.253}},
        {"O2(g)", {154.6, 50.43, 0.022}},
        {"Rn(g)", {377.0, 62.8, 0.0}},
        {"SO2(g)", {430.8, 78.84, 0.245}},
        {"Xe(g)", {289.7, 58.4, 0.0}}};

} // namespace

struct Phreeqc::Impl
{
    // The PHREEQC instance from Phreeqc
    PHREEQC phreeqc;

    // The current temperature in PHREEQC (in units of K)
    double T;

    // The current pressure in PHREEQC (in units of Pa)
    double P;

    // The current ionic strength in PHREEQC (in units of molal)
    double I;

    // The current molar amounts of all species in PHREEQC (in units of mol)
    Vector n;

    // The name of the database file loaded into this instance
    std::string database;

    // The set of elements composing the species
    std::vector<element*> elements;

    // The set of elements composing the species
    std::vector<std::map<element*, double>> elements_in_species;

    // The list of aqueous species (primary and secondary species)
    std::vector<species*> aqueous_species;

    // The list of secondary aqueous species
    std::vector<species*> secondary_species;

    // The set of master aqueous species
    std::vector<species*> master_species;

    // The list of gaseous species
    std::vector<phase*> gaseous_species;

    // The list of mineral species
    std::vector<phase*> mineral_species;

    // The list of reaction equations defining the reactions between product and master species
    std::vector<ReactionEquation> reactions;

    // The index of H2O species in the list of aqueous species
    unsigned iH2O;

    // The elements that are used as redox couples and their redox couples, e.g., N: {N(-3), N(5)}
    std::map<std::string, std::set<std::string>> redox_elements;

    // The names of all elements
    std::vector<std::string> element_names;

    // The names of all species
    std::vector<std::string> species_names;

    // The names of all phases
    std::vector<std::string> phase_names;

    // The formula matrix of the chemical system
    Matrix formula_matrix;

    // The stoichiometric matrix of the equilibrium reactions
    Matrix stoichiometric_matrix;

    // The LU decomposition of the stoichiometric matrix
    Eigen::FullPivLU<Matrix> lu;

    // The electrical charges of the species
    Vector species_charges;

    // The molar masses of the elements (in units of mol/kg)
    Vector element_molar_masses;

    // The ln activity coefficients of the species
    Vector ln_activity_coefficients;

    // The ln activities of the species
    Vector ln_activities;

    // The ln activity constants of the species
    Vector ln_activity_constants;

    // The ln activity coefficients of the aqueous species
    Vector ln_activity_coefficients_aqueous_species;

    // The ln activity coefficients of the gaseous species
    Vector ln_activity_coefficients_gaseous_species;

    // The ln activities of the aqueous species
    Vector ln_activities_aqueous_species;

    // The ln activities of the gaseous species
    Vector ln_activities_gaseous_species;

    // The standard molar Gibbs energies of the species with T and P corrections
    Vector standard_molar_gibbs_energies;

    // The standard molar volumes of the species with T and P corrections
    Vector standard_molar_volumes;

    // The standard molar Gibbs energies of the species with T, P, and I corrections
    Vector standard_molar_gibbs_energies_TPI;

    // The standard molar volumes of the species with T, P, and I corrections
    Vector standard_molar_volumes_TPI;

    // The molar volume of the aqueous phase with T, P, and I corrections
    double molar_volume_aqueous_phase;

    // The molar volume of the gaseous phase with T and P corrections
    double molar_volume_gaseous_phase;

    // The molar volumes of the phases (aqueous, gaseous, minerals, etc.)
    Vector phase_molar_volumes;

    // Construct a default Phreeqc::Impl instance
    Impl();

    // Construct an Phreeqc::Impl instance with given database
    Impl(std::string database);

    // Load a PHREEQC database.
    auto load(std::string database) -> void;

    // Execute a PHREEQC input script file.
    auto execute(std::string input, std::string output) -> void;

    // Execute a PHREEQC input script as a stringstream.
    auto execute(std::stringstream& input, std::string output) -> void;

    // Initialize this Phreeqc::Impl instance according to the active state of PHREEQC
    // This method will check which PHREEQC species and phases are currently active
    // and initialize the data members of this instance with such state.
    auto initialize() -> void;

    // Initialize the species of the chemical system
    auto initializeSpecies() -> void;

    // Initialize the master species that compose the species
    auto initializeMasterSpecies() -> void;

    // Initialize the list of elements that are used as redox couples
    auto initializeRedoxElements() -> void;

    // Initialize the elements that compose the species
    auto initializeElements() -> void;

    // Initialize the names of the elements, species and phases
    auto initializeNames() -> void;

    // Initialize the electrical charges of the species
    auto initializeSpeciesCharges() -> void;

    // Initialize the molar masses of the elements and species
    auto initializeElementMolarMasses() -> void;

    // Initialize the system of reactions between product and master species
    auto initializeReactions() -> void;

    // Initialize the formula matrix of the chemical system
    auto initializeFormulaMatrix() -> void;

    // Initialize the stoichiometric matrix of the chemical system
    auto initializeStoichiometricMatrix() -> void;

    // Initialize the critical properties of the gaseous species
    auto initializeCriticalPropertiesGaseousSpecies() -> void;

    // Initialize the chemical state of the phases and species
    auto initializeChemicalState() -> void;

    // Set the temperature and pressure
    auto set(double T, double P) -> void;

    // Set the temperature, pressure and species composition
    auto set(double T, double P, const Vector& n) -> void;

    // Update the properties with T and P dependency.
    auto updateThermoProperties() -> void;

    // Update the properties with T, P, n dependency.
    auto updateChemicalProperties() -> void;

    // Update the thermodynamic properties of the aqueous phase
    auto updateAqueousProperties() -> void;

    // Update the thermodynamic properties of the gaseous phase
    auto updateGaseousProperties() -> void;

    // Return the number of elements
    auto numElements() const -> unsigned;

    // Return the number of species
    auto numSpecies() const -> unsigned;

    // Return the number of phases
    auto numPhases() const -> unsigned;

    // Return the number of reactions
    auto numReactions() const -> unsigned;

    // Set the molar amounts of the species
    auto setSpeciesAmounts(VectorConstRef n) -> void;

    // Return the temperature of the Phreeqc instance (in units of K)
    auto temperature() const -> double;

    // Return the pressure of the Phreeqc instance (in units of Pa)
    auto pressure() const -> double;

    // Return the molar amounts of the species (in units of mol)
    auto speciesAmounts() const -> VectorConstRef;

    // Return the ln equilibrium constants of the reactions
    auto lnEquilibriumConstants() -> Vector;

    // Return the molar Gibbs energies of the species (in units of J/mol)
    auto speciesMolarGibbsEnergies() -> Vector;

    // Return the molar volumes of the species (in units of J/mol)
    auto speciesMolarVolumes() -> Vector;

    // Return the molar volumes of each phase (in units of m3/mol)
    auto phaseMolarVolumes() -> Vector;

    // Return the natural logarithm of the activity coefficients of the species
    auto lnActivityCoefficients() -> Vector;

    // Return the natural logarithm of the activity constants of the species
    auto lnActivityConstants() -> Vector;

    // Return the natural logarithm of the activities of the species
    auto lnActivities() -> Vector;
};

Phreeqc::Impl::Impl()
{}

Phreeqc::Impl::Impl(std::string database)
    : Impl()
{
    load(database);
}

auto Phreeqc::Impl::load(std::string filename) -> void
{
    //------------------------------------------------------
    // Warning: This method assumes that the Phreeqc::Impl
    // has been reset already.
    //------------------------------------------------------
    // Set the name of the database file
    database = filename;

    // Load the given Phreeqc database
    PhreeqcUtils::load(phreeqc, database);
}

auto Phreeqc::Impl::execute(std::string input, std::string output) -> void
{
    // Execute the given input script file
    PhreeqcUtils::execute(phreeqc, input, output);

    // Initialize the data members after executing the PHREEQC script
    initialize();
}

auto Phreeqc::Impl::initialize() -> void
{
    // Initialize the species pointers
    initializeSpecies();

    // Initialize the corresponding master species of the aqueous species
    initializeMasterSpecies();

    // Initialize the set of elements that are used as redox couples
    initializeRedoxElements();

    // Initialize the set of elements that compose the species
    initializeElements();

    // Initialize the names of the elements, species and phases
    initializeNames();

    // Initialize the electrical charges of the species
    initializeSpeciesCharges();

    // Initialize the molar masses of the elements
    initializeElementMolarMasses();

    // Initialize the system of reactions between product and master species
    initializeReactions();

    // Initialize the formula matrix
    initializeFormulaMatrix();

    // Initialize the stoichiometric matrix
    initializeStoichiometricMatrix();

    // Initialize the critical properties of the gases
    initializeCriticalPropertiesGaseousSpecies();

    // Initialize the chemical state
    initializeChemicalState();
}

auto Phreeqc::Impl::initializeSpecies() -> void
{
    // Initialize the list of all active aqueous species in Phreeqc
    aqueous_species = PhreeqcUtils::activeAqueousSpecies(phreeqc);

    // Initialize the list of secondary aqueous species
    secondary_species = PhreeqcUtils::activeProductSpecies(phreeqc);

    // Initialize the list of gaseous species defined in a gas phase
    gaseous_species = PhreeqcUtils::activeGaseousSpecies(phreeqc);

    // Initialize the list of mineral species active in Phreeqc
    mineral_species = PhreeqcUtils::activePhasesInEquilibriumPhases(phreeqc);

    // Sort the species in alphabetical order
    auto compare_s = [](species* l, species* r) { return std::strcmp(l->name, r->name) < 0; };
    auto compare_p = [](phase* l, phase* r) { return std::strcmp(l->name, r->name) < 0; };

    std::sort(aqueous_species.begin(), aqueous_species.end(), compare_s);
    std::sort(gaseous_species.begin(), gaseous_species.end(), compare_p);
    std::sort(mineral_species.begin(), mineral_species.end(), compare_p);

    // Initialize the index of water among the aqueous species
    iH2O = PhreeqcUtils::index("H2O", aqueous_species);
}

auto Phreeqc::Impl::initializeMasterSpecies() -> void
{
    auto get_corresponding_master_species = [&](species* s) -> species* {
        // Check if the species in e-
        if(s == phreeqc.s_eminus)
            return s;

        for(int i = 0; i < phreeqc.count_species_list; ++i)
            if(phreeqc.species_list[i].s->name == s->name)
                return phreeqc.species_list[i].master_s;

        RuntimeError("Could not initialize the list of master species.",
                     "Could not find the master species of species `" + std::string(s->name) + "`.");

        return nullptr;
    };

    master_species.clear();
    master_species.reserve(phreeqc.count_master);
    for(auto s : aqueous_species)
        master_species.push_back(get_corresponding_master_species(s));
}

auto Phreeqc::Impl::initializeRedoxElements() -> void
{
    redox_elements.clear();
    for(int i = 0; i < phreeqc.count_unknowns; ++i) {
        if(phreeqc.x[i]->master == nullptr)
            continue;
        for(int j = i + 1; j < phreeqc.count_unknowns; ++j) {
            if(phreeqc.x[j]->master == nullptr)
                continue;
            if(phreeqc.x[i]->master[0]->elt->primary == phreeqc.x[j]->master[0]->elt->primary) {
                const auto primary_element = phreeqc.x[i]->master[0]->elt->primary->elt->name;
                const auto secondary_element1 = phreeqc.x[i]->master[0]->elt->name;
                const auto secondary_element2 = phreeqc.x[j]->master[0]->elt->name;

                if(!redox_elements.count(primary_element))
                    redox_elements[primary_element] = {};

                redox_elements[primary_element].insert(secondary_element1);
                redox_elements[primary_element].insert(secondary_element2);
            }
        }
    }
}

auto Phreeqc::Impl::initializeElements() -> void
{
    auto choose_element_is_species = [&](element* e, species* master) {
        if(master == phreeqc.s_eminus)
            return e;
        if(redox_elements.count(e->name))
            if(master->secondary != nullptr &&
               master->secondary->elt->primary->elt->name == std::string(e->name))
                return master->secondary->elt;
            else
                return e->master->s->secondary->elt;
        else
            return e;
    };

    auto choose_element_in_phase = [&](element* e) {
        return redox_elements.count(e->name) ? e->master->s->secondary->elt : e;
    };

    auto get_elements_in_species = [&](species* s, species* master) {
        std::map<element*, double> elements;
        for(auto iter = s->next_elt; iter->elt != nullptr; ++iter)
            elements.insert({choose_element_is_species(iter->elt, master), iter->coef});
        return elements;
    };

    auto get_elements_in_phase = [&](phase* s) {
        std::map<element*, double> elements;
        for(auto iter = s->next_elt; iter->elt != nullptr; ++iter)
            elements.insert({choose_element_in_phase(iter->elt), iter->coef});
        return elements;
    };

    // Clear the member before adding new items to it
    elements_in_species.clear();

    // Collect the elements in the aqueous species
    for(unsigned i = 0; i < aqueous_species.size(); ++i)
        elements_in_species.push_back(get_elements_in_species(aqueous_species[i], master_species[i]));

    // Collect the elements in the gaseous species
    for(auto x : gaseous_species)
        elements_in_species.push_back(get_elements_in_phase(x));

    // Collect the elements in the mineral species
    for(auto x : mineral_species)
        elements_in_species.push_back(get_elements_in_phase(x));

    // Collect all element in a set container
    std::set<element*> element_set;
    for(auto map : elements_in_species)
        for(auto pair : map)
            //            if(pair.first->name != std::string("E") &&
            //                pair.first->name != std::string("e")) // do not add e- as an element
            element_set.insert(pair.first);

    // Transform a std::set to a std::vector of elements
    elements.resize(element_set.size());
    elements.assign(element_set.begin(), element_set.end());

    // Sort the elements in alphabetical order
    std::sort(elements.begin(), elements.end(), [](element* l, element* r) { return std::string(l->name) < std::string(r->name); });
}

auto Phreeqc::Impl::initializeNames() -> void
{
    // Clear container members before adding new items to it
    element_names.clear();
    species_names.clear();
    phase_names.clear();

    // Initialize the names of the elements (alphabetical order)
    for(auto x : elements)
        element_names.push_back(x->name);

    // Add the electrical charge element denoted by Z
    element_names.push_back("Z");

    // Initialize the names of the species (phase order)
    for(auto x : aqueous_species)
        species_names.push_back(x->name);

    for(auto x : gaseous_species)
        species_names.push_back(x->name);

    for(auto x : mineral_species)
        species_names.push_back(x->name);

    // Initialize the names of the phases (phase order)
    phase_names.push_back("Aqueous");

    if(gaseous_species.size())
        phase_names.push_back("Gaseous");

    for(auto x : mineral_species)
        phase_names.push_back(x->name);
}

auto Phreeqc::Impl::initializeSpeciesCharges() -> void
{
    species_charges = zeros(species_names.size());
    unsigned i = 0;
    for(auto x : aqueous_species)
        species_charges[i++] = x->z;
}

auto Phreeqc::Impl::initializeElementMolarMasses() -> void
{
    auto get_element_molar_mass = [](element* e) {
        if(e->name == std::string("E") || e->name == std::string("e"))
            return 0.0; // electron e- element
        return e->primary->elt->gfw;
    };

    const unsigned num_elements = element_names.size();
    element_molar_masses.resize(num_elements);
    for(unsigned i = 0; i < num_elements - 1; ++i) // all, except charge element (last)
        element_molar_masses[i] = get_element_molar_mass(elements[i]) * gram_to_kilogram;
    element_molar_masses[num_elements - 1] = 0.0; // the molar mass of charge element
}

auto Phreeqc::Impl::initializeReactions() -> void
{
    // Clear the reactions member before initializing it
    reactions.clear();

    // Iterate over all aqueous secondary species and get their reaction equation
    for(auto species : secondary_species)
        reactions.push_back(PhreeqcUtils::reactionEquation(species));

    // Iterate over all gaseous species and get their reaction equation
    for(auto species : gaseous_species)
        reactions.push_back(PhreeqcUtils::reactionEquation(species));

    // Iterate over all pure mineral species and get their reaction equation
    for(auto species : mineral_species)
        reactions.push_back(PhreeqcUtils::reactionEquation(species));
}

auto Phreeqc::Impl::initializeFormulaMatrix() -> void
{
    const unsigned num_elements = element_names.size();
    const unsigned num_species = species_names.size();

    auto get_element_stoichiometry = [&](unsigned ielement, unsigned ispecies) {
        // Check if the element index corresponds to the charge index (last index)
        if(ielement == element_names.size() - 1)
            return species_charges[ispecies];
        auto element = elements[ielement];
        for(auto pair : elements_in_species[ispecies])
            if(pair.first == element)
                return pair.second;
        return 0.0;
    };

    formula_matrix.resize(num_elements, num_species);
    for(unsigned i = 0; i < num_species; ++i)
        for(unsigned j = 0; j < num_elements; ++j)
            formula_matrix(j, i) = get_element_stoichiometry(j, i);
}

auto Phreeqc::Impl::initializeStoichiometricMatrix() -> void
{
    // Define the number of reactions and species
    const unsigned num_reactions = reactions.size();
    const unsigned num_species = numSpecies();

    // Initialize the stoichiometric matrix of the equilibrium reactions
    stoichiometric_matrix = zeros(num_reactions, num_species);
    for(unsigned j = 0; j < num_reactions; ++j) {
        for(auto pair : reactions[j]) {
            const std::string species_name = pair.first;
            const double species_coef = pair.second;
            const unsigned i = index(species_name, species_names);
            stoichiometric_matrix(j, i) = species_coef;
        }
    }

    // Initialize the LU decomposition of the stoichiometric matrix
    lu.compute(stoichiometric_matrix);
}

auto Phreeqc::Impl::initializeCriticalPropertiesGaseousSpecies() -> void
{
    // Iterate over the gaseous species and set their critical properties.
    // However, only do this for those gases whose critical properties were
    // not specified in the PHASES block in the PHREEQC script file.
    for(auto species : gaseous_species) {
        auto iter = critical_properties.find(species->name);

        if(iter != critical_properties.end()) {
            if(species->t_c == 0)
                species->t_c = iter->second[0];
            if(species->p_c == 0)
                species->p_c = convertBarToAtm(iter->second[1]);
            if(species->omega == 0)
                species->omega = iter->second[2];
        }
    }
}

void Phreeqc::Impl::initializeChemicalState()
{
    // Initialize temperature and pressure
    T = phreeqc.tk_x;
    P = phreeqc.patm_x * atm_to_pascal;

    // Initialize the amounts of species
    n.resize(numSpecies());

    n << PhreeqcUtils::speciesAmounts(phreeqc, aqueous_species),
        PhreeqcUtils::speciesAmounts(phreeqc, gaseous_species),
        PhreeqcUtils::speciesAmounts(phreeqc, mineral_species);

    // Initialize thermodynamic and chemical properties
    set(T, P, n);
}

auto Phreeqc::Impl::set(double T, double P) -> void
{
    // Set the current temperature and pressure of PHREEQC
    this->T = T;
    this->P = P;

    // Update the temperature member (in units of K)
    phreeqc.tk_x = T;

    // Update the temperature member (in units of celsius)
    phreeqc.tc_x = T - 273.15;

    // Update the pressure member (in units of atm)
    phreeqc.patm_x = P * pascal_to_atm;

    // Update the thermodynamic properties (T and P dependent)
    updateThermoProperties();
}

auto Phreeqc::Impl::set(double T, double P, const Vector& n) -> void
{
    // Set the current temperature and pressure of PHREEQC
    this->T = T;
    this->P = P;

    // Update the temperature member (in units of K)
    phreeqc.tk_x = T;

    // Update the temperature member (in units of celsius)
    phreeqc.tc_x = T - 273.15;

    // Update the pressure member (in units of atm)
    phreeqc.patm_x = P * pascal_to_atm;

    // Update the amounts of the species
    setSpeciesAmounts(n);

    // Update the thermodynamic properties (T and P dependent)
    updateThermoProperties();

    // Update the thermodynamic properties (T, P, and n dependent)
    updateChemicalProperties();
}

auto Phreeqc::Impl::updateThermoProperties() -> void
{
    // Set ionic strength to zero to eliminate ionic strength corrections in
    // the equilibrium constants of reactions. These ionic strength corrections
    // are used as contributions in the activities of the species, which are
    // then evaluated whenever composition changes.
    phreeqc.mu_x = 0.0;

    // Update equilibrium constants of reactions with T and P corrections.
    // This Phreeqc::k_temp call also updates density and dielectric
    // properties of water at the given T and P conditions.
    phreeqc.k_temp(phreeqc.tc_x, phreeqc.patm_x);

    // Set the ionic strength back to its original value
    phreeqc.mu_x = I;

    // Update the standard Gibbs energies of the species with T and P corrections
    standard_molar_gibbs_energies = speciesMolarGibbsEnergies();

    // Update the standard molar volumes of the species with T and P corrections
    standard_molar_volumes = speciesMolarVolumes();

    // Update the ln activity constants
    ln_activity_constants = lnActivityConstants();
}

auto Phreeqc::Impl::updateChemicalProperties() -> void
{
    // Update equilibrium constants of reactions with T, P, and I corrections.
    // This Phreeqc::k_temp call also updates density and dielectric
    // properties of water at the given T and P conditions.
    phreeqc.k_temp(phreeqc.tc_x, phreeqc.patm_x);

    // Update the standard Gibbs energies of the species with T, P, and I corrections
    standard_molar_gibbs_energies_TPI = speciesMolarGibbsEnergies();

    // Update the standard molar volumes of the species with T, P, and I corrections
    standard_molar_volumes_TPI = speciesMolarVolumes();

    updateAqueousProperties();

    updateGaseousProperties();

    ln_activity_coefficients = lnActivityCoefficients();

    ln_activities = lnActivities();

    phase_molar_volumes = phaseMolarVolumes();
}

auto Phreeqc::Impl::updateAqueousProperties() -> void
{
    // Define some auxiliary variables
    const unsigned num_aqueous_species = aqueous_species.size();
    const double ln_10 = std::log(10.0);

    // Define some auxiliary alias
    Vector& ln_g = ln_activity_coefficients_aqueous_species;
    Vector& ln_a = ln_activities_aqueous_species;

    // Allocate memory
    ln_g.resize(num_aqueous_species);
    ln_a.resize(num_aqueous_species);

    // Calculate the activity coefficients of the aqueous species
    if(phreeqc.pitzer_model || phreeqc.sit_model) {
        // Calculate the activity coefficients using either Pitzer or SIT models
        if(phreeqc.pitzer_model)
            phreeqc.pitzer();
        else
            phreeqc.sit();

        // Collect the updated activity coefficients
        unsigned ispecies = 0;
        for(auto species : aqueous_species)
            ln_g[ispecies++] = species->lg_pitzer * ln_10;
    } else {
        // Calculate the activity coefficients using conventional Phreeqc models
        phreeqc.gammas(phreeqc.mu_x);

        // Collect the updated activity coefficients
        unsigned ispecies = 0;
        for(auto species : aqueous_species)
            ln_g[ispecies++] = species->lg * ln_10;
    }

    // Calculate the natural log of the activities of the aqueous species
    for(unsigned i = 0; i < num_aqueous_species; ++i)
        ln_a[i] = std::isfinite(aqueous_species[i]->lm) ? ln_g[i] + aqueous_species[i]->lm * ln_10 : 0.0;

    // Get the molar amounts of the aqueous species
    const auto n_aqueous = n.head(num_aqueous_species);

    // Calculate the total amount of moles in the aqueous phase
    const double n_total = sum(n_aqueous);

    // Calculate the mole fraction of H2O
    const double nH2O = aqueous_species[iH2O]->moles;
    const double xH2O = nH2O / n_total;

    // Calculate the activity of water
    if(phreeqc.pitzer_model || phreeqc.sit_model)
        ln_a[iH2O] = std::log(phreeqc.AW);
    else
        ln_a[iH2O] = std::log(xH2O);

    // Get the molar volumes of the aqueous species with T, P, I corrections
    const auto v_aqueous = standard_molar_volumes_TPI.head(num_aqueous_species);

    // Calculate the molar volume of the aqueous phase with T, P, I corrections
    molar_volume_aqueous_phase = (n_total > 0) ? dot(v_aqueous, n_aqueous) / n_total : 0.0;
}

auto Phreeqc::Impl::updateGaseousProperties() -> void
{
    // The number of aqueous and gaseous species
    const unsigned num_aqueous_species = aqueous_species.size();
    const unsigned num_gaseous_species = gaseous_species.size();

    // Skip this function if the system has no gaseous phase
    if(num_gaseous_species == 0)
        return;

    // Define some auxiliary variables
    const double T = temperature();
    const double P = pressure();
    const double Patm = P * pascal_to_atm;
    const double Pbar = P * pascal_to_bar;

    // Define some auxiliary alias
    auto& ln_g = ln_activity_coefficients_gaseous_species;
    auto& ln_a = ln_activities_gaseous_species;

    // Allocate memory
    ln_g.resize(num_gaseous_species);
    ln_a.resize(num_gaseous_species);

    // Calculate the thermodynamic properties of the gaseous phase using Peng-Robinson EOS
    const double v = phreeqc.calc_PR(gaseous_species, Patm, T, 0.0);

    // Calculate the ln activity coefficients and ln activities of the gaseous species
    for(unsigned i = 0; i < num_gaseous_species; ++i) {
        const double x = gaseous_species[i]->fraction_x; // the molar fraction of the gas
        const double phi = gaseous_species[i]->pr_phi;   // the fugacity coefficient of the gas
        ln_g[i] = std::log(phi);
        ln_a[i] = std::log(x * phi * Pbar);
    }

    // Get the molar amounts of the gaseous species
    const auto n_gaseous = n.segment(num_aqueous_species, num_gaseous_species);

    // Calculate the total amount of moles in the gaseous phase
    const double n_total = sum(n_gaseous);

    // Ensure the molar volume of the phase is zero if it has zero moles
    molar_volume_gaseous_phase = (n_total > 0.0) ? v : 0.0;
}

auto Phreeqc::Impl::numElements() const -> unsigned
{
    return element_names.size();
}

auto Phreeqc::Impl::numSpecies() const -> unsigned
{
    return species_names.size();
}

auto Phreeqc::Impl::numPhases() const -> unsigned
{
    return phase_names.size();
}

auto Phreeqc::Impl::numReactions() const -> unsigned
{
    return secondary_species.size() + gaseous_species.size() + mineral_species.size();
}

auto Phreeqc::Impl::setSpeciesAmounts(VectorConstRef n) -> void
{
    // Get the number of aqueous, gaseous and mineral species
    const unsigned num_aqueous = aqueous_species.size();
    const unsigned num_gaseous = gaseous_species.size();
    const unsigned num_mineral = mineral_species.size();

    // Get the segments of aqueous, gaseous and mineral species
    auto n_aqueous = n.head(num_aqueous);
    auto n_gaseous = n.segment(num_aqueous, num_gaseous);
    auto n_mineral = n.tail(num_mineral);

    // Get data related to water
    const double nH2O = n_aqueous[iH2O];
    const double massH2O = nH2O * waterMolarMass;

    // Set the copy of current molar amounts of all species in PHREEQC
    this->n = n;

    // Set the molar amounts of species and phase instances of PHREEQC
    // for calculation of phase properties such as activity coefficients,
    // phase densities, etc.

    // Set the molar amounts of the aqueous species
    for(unsigned i = 0; i < num_aqueous; ++i)
        aqueous_species[i]->moles = n_aqueous[i];

    // Set the molar amounts of the gaseous species
    for(unsigned i = 0; i < num_gaseous; ++i)
        gaseous_species[i]->moles_x = n_gaseous[i];

    // Set the molar amounts of the mineral species
    for(unsigned i = 0; i < num_mineral; ++i)
        mineral_species[i]->moles_x = n_mineral[i];

    // Check if mass of water is positive before updating aqueous properties
    if(massH2O > 0.0) {
        // Set the molalities of the aqueous species
        for(unsigned i = 0; i < num_aqueous; ++i)
            aqueous_species[i]->lm = std::log10(n_aqueous[i] / massH2O);

        // Update the ionic strength of the aqueous phase
        I = 0.0;
        for(auto species : aqueous_species)
            I += species->moles * species->z * species->z;
        I *= 0.5 / massH2O;

        // Set the ionic strength of the aqueous solution
        phreeqc.mu_x = I;
    } else {
        // Ionic strength is zero as a result of zero amounts of water
        I = phreeqc.mu_x = 0.0;
    }
}

auto Phreeqc::Impl::temperature() const -> double
{
    return T;
}

auto Phreeqc::Impl::pressure() const -> double
{
    return P;
}

auto Phreeqc::Impl::speciesAmounts() const -> VectorConstRef
{
    return n;
}

auto Phreeqc::Impl::lnEquilibriumConstants() -> Vector
{
    // The natural log of 10
    const auto ln10 = 2.30258509299;

    // Collect the ln equilibrium constants of the reactions
    Vector ln_k(numReactions());

    auto ireaction = 0;
    for(auto species : secondary_species)
        ln_k[ireaction++] = species->lk * ln10;

    for(auto species : gaseous_species)
        ln_k[ireaction++] = species->lk * ln10;

    for(auto species : mineral_species)
        ln_k[ireaction++] = species->lk * ln10;

    return ln_k;
}

auto Phreeqc::Impl::speciesMolarVolumes() -> Vector
{
    // Define some auxiliary variables
    const auto R = universalGasConstant;

    // Collect the standard molar volumes of the species
    Vector v(numSpecies());

    auto ispecies = 0;
    for(auto species : aqueous_species)
        v[ispecies++] = species->logk[vm_tc] * cm3_to_m3;

    for(auto species : gaseous_species)
        v[ispecies++] = R * T / P;

    for(auto species : mineral_species)
        v[ispecies++] = species->logk[vm0] * cm3_to_m3;

    return v;
}

auto Phreeqc::Impl::speciesMolarGibbsEnergies() -> Vector
{
    // The universal gas constant (in units of J/(mol*K))
    const double R = universalGasConstant;

    // Calculate the natural log of the equilibrium constants
    Vector ln_k = lnEquilibriumConstants();

    // Use the LU decomposition of the stoichiometric matrix to calculate `u0`
    Vector G = lu.solve(ln_k);
    G *= -R * T;

    return G;
}

auto Phreeqc::Impl::phaseMolarVolumes() -> Vector
{
    const unsigned num_phases = numPhases();

    Vector vphases(num_phases);

    if(num_phases > 0) {
        unsigned offset = 0;

        if(aqueous_species.size())
            vphases[offset++] = molar_volume_aqueous_phase;

        if(gaseous_species.size())
            vphases[offset++] = molar_volume_gaseous_phase;

        vphases.tail(num_phases - offset) = standard_molar_volumes.tail(num_phases - offset);
    }

    return vphases;
}

auto Phreeqc::Impl::lnActivityCoefficients() -> Vector
{
    const unsigned num_minerals = mineral_species.size();

    Vector ln_g(numSpecies());
    ln_g << ln_activity_coefficients_aqueous_species,
        ln_activity_coefficients_gaseous_species,
        zeros(num_minerals);

    return ln_g;
}

auto Phreeqc::Impl::lnActivityConstants() -> Vector
{
    // The number of aqueous, gaseous, and mineral species
    const unsigned num_aqueous_species = aqueous_species.size();
    const unsigned num_gaseous_species = gaseous_species.size();
    const unsigned num_mineral_species = mineral_species.size();

    // Convert pressure from pascal to bar
    const double Pbar = 1e-5 * P;

    Vector res(numSpecies());
    res << constants(num_aqueous_species, std::log(55.508472)),
        constants(num_gaseous_species, std::log(Pbar)),
        zeros(num_mineral_species);

    // The ln activity constant for water is zero
    res[iH2O] = 0.0;

    return res;
}

auto Phreeqc::Impl::lnActivities() -> Vector
{
    // Auxiliary variables
    const auto R = universalGasConstant;
    const auto num_minerals = mineral_species.size();

    // Collect the activities of the species
    Vector ln_a(numSpecies());
    ln_a << ln_activities_aqueous_species,
        ln_activities_gaseous_species,
        zeros(num_minerals);

    // The standard molar Gibbs energies with T and P corrections
    const auto& G0 = standard_molar_gibbs_energies;

    // The standard molar Gibbs energies with T, P, and I corrections
    const auto& G0TPI = standard_molar_gibbs_energies_TPI;

    // Here is where the ionic strength corrections on the equilibrium constants
    // of reactions are transferred to the activities of the species.
    // This is needed because standard properties should not depend on composition,
    // but only on T and P.
    ln_a += (G0TPI - G0) * 1.0 / (R * T);

    return ln_a;
}

Phreeqc::Phreeqc()
    : pimpl(new Impl())
{}

Phreeqc::Phreeqc(std::string database)
    : pimpl(new Impl(database))
{}

Phreeqc::~Phreeqc()
{}

auto Phreeqc::temperature() const -> double
{
    return pimpl->temperature();
}

auto Phreeqc::pressure() const -> double
{
    return pimpl->pressure();
}

auto Phreeqc::speciesAmounts() const -> Vector
{
    return pimpl->speciesAmounts();
}

auto Phreeqc::numElements() const -> unsigned
{
    return pimpl->numElements();
}

auto Phreeqc::numSpecies() const -> unsigned
{
    return pimpl->numSpecies();
}

auto Phreeqc::numPhases() const -> unsigned
{
    return pimpl->numPhases();
}

auto Phreeqc::numSpeciesInPhase(Index iphase) const -> unsigned
{
    const unsigned num_aqueous = pimpl->aqueous_species.size();
    const unsigned num_gaseous = pimpl->gaseous_species.size();

    // Ensure `index` is not out-of-bound
    Assert(iphase < numPhases(),
           "Cannot get the number of species in phase with index "
           "`" +
               std::to_string(iphase) + "`.",
           "The given index is out of range.")

        // Return the number of aqueous species
        if(iphase == 0) return num_aqueous;

    // Return the number of gaseous species if they exist
    if(iphase == 1 && num_gaseous)
        return num_gaseous;

    // Return 1 as all mineral phases only have one species
    return 1;
}

auto Phreeqc::elementName(Index ielement) const -> std::string
{
    return pimpl->element_names[ielement];
}

auto Phreeqc::elementMolarMass(Index ielement) const -> double
{
    return pimpl->element_molar_masses[ielement];
}

auto Phreeqc::elementStoichiometry(Index ispecies, Index ielement) const -> double
{
    return pimpl->formula_matrix(ielement, ispecies);
}

auto Phreeqc::speciesName(Index ispecies) const -> std::string
{
    return pimpl->species_names[ispecies];
}

auto Phreeqc::phaseName(Index iphase) const -> std::string
{
    return pimpl->phase_names[iphase];
}

auto Phreeqc::set(double T, double P) -> void
{
    pimpl->set(T, P);
}

auto Phreeqc::set(double T, double P, VectorConstRef n) -> void
{
    pimpl->set(T, P, n);
}

auto Phreeqc::load(std::string database) -> void
{
    // Resets this instance before loading a new database file
    reset();
    pimpl->load(database);
}

auto Phreeqc::execute(std::string input, std::string output) -> void
{
    pimpl->execute(input, output);
}

auto Phreeqc::execute(std::string input) -> void
{
    pimpl->execute(input, {});
}

auto Phreeqc::reset() -> void
{
    pimpl.reset(new Phreeqc::Impl());
}

auto Phreeqc::reactions() const -> const std::vector<ReactionEquation>&
{
    return pimpl->reactions;
}

auto Phreeqc::stoichiometricMatrix() const -> Matrix
{
    return pimpl->stoichiometric_matrix;
}

auto Phreeqc::standardMolarGibbsEnergies() const -> Vector
{
    return pimpl->standard_molar_gibbs_energies;
}

auto Phreeqc::standardMolarEnthalpies() const -> Vector
{
    return zeros(numSpecies());
}

auto Phreeqc::standardMolarVolumes() const -> Vector
{
    return pimpl->standard_molar_volumes;
}

auto Phreeqc::standardMolarHeatCapacitiesConstP() const -> Vector
{
    return zeros(numSpecies());
}

auto Phreeqc::standardMolarHeatCapacitiesConstV() const -> Vector
{
    return zeros(numSpecies());
}

auto Phreeqc::lnActivityCoefficients() const -> Vector
{
    return pimpl->lnActivityCoefficients();
}

auto Phreeqc::lnActivityConstants() const -> Vector
{
    return pimpl->lnActivityConstants();
}

auto Phreeqc::lnActivities() const -> Vector
{
    return pimpl->lnActivities();
}

auto Phreeqc::lnEquilibriumConstants() const -> Vector
{
    return pimpl->lnEquilibriumConstants();
}

auto Phreeqc::phaseMolarVolumes() const -> Vector
{
    return pimpl->phaseMolarVolumes();
}

auto Phreeqc::properties(ThermoModelResult& res, double T, double P) -> void
{
    // Update the temperature and pressure of the Phreeqc instance
    set(T, P);

    // Set the thermodynamic properties of the phases and species
    res.standardPartialMolarGibbsEnergies().val = pimpl->standard_molar_gibbs_energies;
    res.standardPartialMolarVolumes().val = pimpl->standard_molar_volumes;
    res.lnActivityConstants().val = pimpl->ln_activity_constants;
}

auto Phreeqc::properties(ChemicalModelResult& res, double T, double P, VectorConstRef n) -> void
{
    // Update the temperature, pressure, and species amounts of the Phreeqc instance
    set(T, P, n);

    // Set the molar volumes of the phases, ln activity coefficients and ln activities of all species
    res.phaseMolarVolumes().val = pimpl->phase_molar_volumes;
    res.lnActivityCoefficients().val = pimpl->ln_activity_coefficients;
    res.lnActivities().val = pimpl->ln_activities;

    // The number of phases
    const Index num_phases = numPhases();

    // Set the molar derivatives of the activities
    Index offset = 0;
    for(Index iphase = 0; iphase < num_phases; ++iphase) {
        // The number of species in the current phase
        const Index size = numSpeciesInPhase(iphase);

        // The species amounts in the current phase
        const auto np = n.segment(offset, size);

        // Set d(ln(a))/dn to d(ln(x))/dn, where x is mole fractions
        res.lnActivities().ddn.block(offset, offset, size, size) = -1.0 / sum(np) * ones(size, size);
        res.lnActivities().ddn.block(offset, offset, size, size).diagonal() += 1.0 / np;

        offset += size;
    }
}

auto Phreeqc::clone() const -> std::shared_ptr<Interface>
{
    return std::make_shared<Phreeqc>(*this);
}

auto Phreeqc::phreeqc() -> PHREEQC&
{
    return pimpl->phreeqc;
}

auto Phreeqc::phreeqc() const -> const PHREEQC&
{
    return pimpl->phreeqc;
}

} // namespace Reaktoro
