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

#include "Phreeqx.hpp"

// Eigen includes
#include <eigen/Dense>

// Reaktoro includes
#include "internal/PhreeqcUtils.hpp"
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

// Phreeqc includes
#define protected public
#include <phreeqc/Phreeqc.h>
#include <phreeqc/GasPhase.h>

namespace Reaktoro {
namespace {

const double gram_to_kilogram = 1e-3;
const double kilogram_to_gram = 1e+3;

const double pascal_to_bar = 1e-5;
const double bar_to_pascal = 1e+5;

const double pascal_to_atm = 9.86923267e-6;
const double atm_to_pascal = 1/pascal_to_atm;

const double m3_to_liter = 1e+3;
const double liter_to_m3 = 1e-3;

const double m3_to_cm3 = 1e+6;
const double cm3_to_m3 = 1e-6;

} // namespace

struct Phreeqx::Impl
{
    /// The TNode instance from Phreeqc
    Phreeqc phreeqc;

    // The set of elements composing the species
    std::vector<element*> elements;

    // The list of aqueous species (primary and secondary species)
    std::vector<species*> aqueous_species;

    // The list of secondary aqueous species
    std::vector<species*> secondary_species;

    // The list of gaseous species
    std::vector<phase*> gaseous_species;

    // The list of mineral species
    std::vector<phase*> mineral_species;

    // The index of H2O species in the list of aqueous species
    unsigned iH2O;

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

    // The SVD decomposition of the stoichiometric matrix
    Eigen::JacobiSVD<Matrix> svd;

    // The molar masses of the elements (in units of mol/kg)
    Vector element_molar_masses;

    // Construct a default Impl instance
    Impl();

    // Construct a custom Impl instance
    Impl(std::string database, std::string script);

    // Initialize the species of the chemical system
    auto initializeSpecies() -> void;

    // Initialize the elements that compose the species
    auto initializeElements() -> void;

    // Initialize the names of the elements, species and phases
    auto initializeNames() -> void;

    // Initialize the formula matrix of the chemical system
    auto initializeFormulaMatrix() -> void;

    // Initialize the stoichiometric matrix of the chemical system
    auto initializeStoichiometricMatrix() -> void;

    // Initialize the molar masses of the elements and species
    auto initializeMolarMasses() -> void;

    // Return the number of elements
    auto numElements() const -> unsigned;

    // Return the number of species
    auto numSpecies() const -> unsigned;

    // Return the number of phases
    auto numPhases() const -> unsigned;

    // Return the number of reactions
    auto numReactions() const -> unsigned;

    // Set the molar amounts of the species
    auto setSpeciesAmounts(const Vector& n) -> void;

    /// Return the temperature of the Phreeqc instance (in units of K)
    auto temperature() const -> double;

    /// Return the pressure of the Phreeqc instance (in units of Pa)
    auto pressure() const -> double;

    /// Return the molar amounts of the aqueous species (in units of mol)
    auto speciesAmountsAqueousSpecies() const -> Vector;

    /// Return the molar amounts of the gaseous species (in units of mol)
    auto speciesAmountsGaseousSpecies() const -> Vector;

    /// Return the molar amounts of the mineral species (in units of mol)
    auto speciesAmountsMineralSpecies() const -> Vector;

    /// Return the molar amounts of the species (in units of mol)
    auto speciesAmounts() const -> Vector;

    /// Return the molar amount of a species (in units of mol)
    auto speciesAmount(unsigned index) const -> double;

    // Return the natural logarithm of the equilibrium constants of the reactions
    auto lnEquilibriumConstants() -> Vector;

    // Return the natural logarithm of the activity coefficients of the aqueous species
    auto lnActivityCoefficientsAqueousSpecies() -> Vector;

    // Return the natural logarithm of the activities of the aqueous species
    auto lnActivitiesAqueousSpecies() -> Vector;

    // Return the natural logarithm of the activities of the gaseous species
    auto lnActivitiesGaseousSpecies() -> Vector;

    // Return the natural logarithm of the activities of the mineral species
    auto lnActivitiesMineralSpecies() -> Vector;

    // Return the natural logarithm of the activities of the species
    auto lnActivities() -> Vector;

    // Return the standard molar Gibbs energies of the species (in units of J/mol)
    auto standardMolarGibbsEnergies() -> Vector;

    // Return the standard molar volumes of the aqueous species (in units of m3/mol)
    auto standardMolarVolumesAqueousSpecies() -> Vector;

    // Return the standard molar volumes of the gaseous species (in units of m3/mol)
    auto standardMolarVolumesGaseousSpecies() -> Vector;

    // Return the standard molar volumes of the mineral species (in units of m3/mol)
    auto standardMolarVolumesMineralSpecies() -> Vector;

    // Return the standard molar volumes of the species (in units of m3/mol)
    auto standardMolarVolumes() -> Vector;

    /// Return the molar volume of the aqueous phase (in units of m3/mol)
    auto molarVolumeAqueousPhase() -> double;

    /// Return the molar volume of the gaseous phase (in units of m3/mol)
    auto molarVolumeGaseousPhase() -> double;

    /// Return the molar volumes of the mineral phases (in units of m3/mol)
    auto molarVolumeMineralPhases() -> Vector;

    /// Return the molar volumes of each phase (in units of m3/mol)
    auto phaseMolarVolumes() -> Vector;
};

Phreeqx::Impl::Impl()
{}

Phreeqx::Impl::Impl(std::string database, std::string script)
{
    // Initialize the low-level Phreeqc instance
    loadDatabase(phreeqc, database);
    loadScript(phreeqc, script);

    // Initialize the species pointers
    initializeSpecies();

    // Initialize the set of elements that compose the species
    initializeElements();

    // Initialize the names of the elements, species and phases
    initializeNames();

    // Initialize the formula matrix
    initializeFormulaMatrix();

    // Initialize the stoichiometric matrix
    initializeStoichiometricMatrix();

    // Initialize the molar masses of the elements and species
    initializeMolarMasses();
}

auto Phreeqx::Impl::initializeSpecies() -> void
{
    // Initialize the list of all active aqueous species in Phreeqc
    aqueous_species = collectAqueousSpecies(phreeqc);

    // Initialize the list of secondary aqueous species
    secondary_species = collectSecondarySpecies(phreeqc);

    // Initialize the list of gaseous species defined in a gas phase
    gaseous_species = collectGaseousSpecies(phreeqc);

    // Initialize the list of mineral species active in Phreeqc
    mineral_species = collectMineralSpecies(phreeqc);

    // Initialize the index of water among the aqueous species
    iH2O = index("H2O", aqueous_species);
}

auto Phreeqx::Impl::initializeElements() -> void
{
    std::set<element*> element_set;

    // Collect the elements in the aqueous species
    for(auto x : aqueous_species)
        for(auto y : getElementsInSpecies(x))
            element_set.insert(y.first);

    // Collect the elements in the gaseous species
    for(auto x : gaseous_species)
        for(auto y : getElementsInPhase(x))
            element_set.insert(y.first);

    // Collect the elements in the mineral species
    for(auto x : mineral_species)
        for(auto y : getElementsInPhase(x))
            element_set.insert(y.first);

    // Transform a std::set to a std::vector of elements
    elements.resize(element_set.size());
    elements.assign(element_set.begin(), element_set.end());

    // Sort the elements in alphabetical order
    std::sort(elements.begin(), elements.end(),
        [](element* l, element* r) { return std::strcmp(l->name, r->name); });
}

auto Phreeqx::Impl::initializeNames() -> void
{
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

auto Phreeqx::Impl::initializeFormulaMatrix() -> void
{
    const unsigned num_elements = element_names.size();
    const unsigned num_species = species_names.size();

    formula_matrix.resize(num_elements, num_species);

    unsigned ispecies = 0;
    for(auto species : aqueous_species)
    {
        for(unsigned j = 0; j < num_elements; ++j)
            formula_matrix(j, ispecies) =
                elementStoichiometryInSpecies(element_names[j], species);
        ++ispecies;
    }

    for(auto species : gaseous_species)
    {
        for(unsigned j = 0; j < num_elements; ++j)
            formula_matrix(j, ispecies) =
                elementStoichiometryInPhase(element_names[j], species);
        ++ispecies;
    }

    for(auto species : mineral_species)
    {
        for(unsigned j = 0; j < num_elements; ++j)
            formula_matrix(j, ispecies) =
                elementStoichiometryInPhase(element_names[j], species);
        ++ispecies;
    }
}

auto Phreeqx::Impl::initializeStoichiometricMatrix() -> void
{
    std::vector<std::map<std::string, double>> equations;

    // Iterate over all aqueous secondary species and get their reaction equation
    for(auto species : secondary_species)
        equations.push_back(getReactionEquation(species));

    // Iterate over all gaseous species and get their reaction equation
    for(auto species : gaseous_species)
        equations.push_back(getReactionEquation(species));

    // Iterate over all pure mineral species and get their reaction equation
    for(auto species : mineral_species)
        equations.push_back(getReactionEquation(species));

    // Define the number of reactions and species
    const unsigned num_reactions = equations.size();
    const unsigned num_species = numSpecies();

    // Initialize the stoichiometric matrix of the equilibrium reactions
    stoichiometric_matrix = Matrix::Zero(num_reactions, num_species);
    for(unsigned j = 0; j < num_reactions; ++j)
    {
        for(auto pair : equations[j])
        {
            const std::string species_name = pair.first;
            const double species_coef = pair.second;
            const unsigned i = index(species_name, species_names);
            stoichiometric_matrix(j, i) = species_coef;
        }
    }

    // Initialize the SVD decomposition of the stoichiometric matrix
    svd.compute(stoichiometric_matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
}

auto Phreeqx::Impl::initializeMolarMasses() -> void
{
    const unsigned num_elements = element_names.size();
    element_molar_masses.resize(num_elements);
    for(unsigned i = 0; i < num_elements - 1; ++i) // all, except charge element (last)
        element_molar_masses[i] = elements[i]->gfw * gram_to_kilogram;
    element_molar_masses[num_elements - 1] = 0.0; // the molar mass of charge element
}

auto Phreeqx::Impl::numElements() const -> unsigned
{
    return element_names.size();
}

auto Phreeqx::Impl::numSpecies() const -> unsigned
{
    return species_names.size();
}

auto Phreeqx::Impl::numPhases() const -> unsigned
{
    return phase_names.size();
}

auto Phreeqx::Impl::numReactions() const -> unsigned
{
    return secondary_species.size() + gaseous_species.size() + mineral_species.size();
}

auto Phreeqx::Impl::setSpeciesAmounts(const Vector& n) -> void
{
    // Get the number of aqueous, gaseous and mineral species
    const unsigned num_aqueous = aqueous_species.size();
    const unsigned num_gaseous = gaseous_species.size();
    const unsigned num_mineral = mineral_species.size();

    // Get the segments of aqueous, gaseous and mineral species
    auto n_aqueous = n.topRows(num_aqueous);
    auto n_gaseous = n.middleRows(num_aqueous, num_gaseous);
    auto n_mineral = n.bottomRows(num_mineral);

    // Get data related to water
    const double nH2O = n_aqueous[iH2O];
    const double massH2O = nH2O * waterMolarMass;

    // Set the molar amounts and molalities of the aqueous species
    for(unsigned i = 0; i < num_aqueous; ++i)
    {
        aqueous_species[i]->moles = n_aqueous[i];
        aqueous_species[i]->lm = std::log10(n_aqueous[i]/massH2O);
    }

    // Set the molar amounts of the gaseous species
    for(unsigned i = 0; i < num_gaseous; ++i)
        gaseous_species[i]->moles_x = n_gaseous[i];

    // Set the molar amounts of the mineral species
    for(unsigned i = 0; i < num_mineral; ++i)
        mineral_species[i]->moles_x = n_mineral[i];

    // Calculate the ionic strength of the aqueous phase
    double ionic_strength = 0.0;
    for(auto species : aqueous_species)
        ionic_strength += species->moles * species->z * species->z;
    ionic_strength *= 0.5/massH2O;

    // Set the ionic strength of the aqueous solution
    phreeqc.mu_x = ionic_strength;
}

auto Phreeqx::Impl::temperature() const -> double
{
    return phreeqc.tk_x;
}

auto Phreeqx::Impl::pressure() const -> double
{
    return phreeqc.patm_x * atm_to_pascal;
}

auto Phreeqx::Impl::speciesAmountsAqueousSpecies() const -> Vector
{
    return speciesAmountsInSpecies(aqueous_species);
}

auto Phreeqx::Impl::speciesAmountsGaseousSpecies() const -> Vector
{
    return speciesAmountsInPhases(gaseous_species);
}

auto Phreeqx::Impl::speciesAmountsMineralSpecies() const -> Vector
{
    return speciesAmountsInPhases(mineral_species);
}

auto Phreeqx::Impl::speciesAmounts() const -> Vector
{
    Vector n(numSpecies());
    n << speciesAmountsAqueousSpecies(),
         speciesAmountsGaseousSpecies(),
         speciesAmountsMineralSpecies();
    return n;
}

auto Phreeqx::Impl::speciesAmount(unsigned index) const -> double
{
    const unsigned num_aqueous = aqueous_species.size();
    const unsigned num_gaseous = gaseous_species.size();
    const unsigned num_mineral = mineral_species.size();

    Assert(index < num_aqueous + num_gaseous + num_mineral,
        "Cannot get the amount of species with index `" + std::to_string(index) + "`.",
        "The given index is out of range.");

    // Check if `index` points to an aqueous species
    if(index < num_aqueous)
        return aqueous_species[index]->moles;

    // Check if `index` points to a gaseous species
    if(index < num_aqueous + num_gaseous)
        return gaseous_species[index-num_aqueous]->moles_x;

    // Then `index` must point to a mineral species
    return mineral_species[index-num_aqueous-num_gaseous]->moles_x;
}

auto Phreeqx::Impl::lnEquilibriumConstants() -> Vector
{
    const unsigned num_reactions = numReactions();

    const double T = temperature();
    const double P = pressure();

    Vector ln_k(num_reactions);

    unsigned ireaction = 0;
    for(auto species : secondary_species)
        ln_k[ireaction++] = -lnEquilibriumConstant(species, T, P);

    for(auto species : gaseous_species)
        ln_k[ireaction++] = lnEquilibriumConstant(species, T, P);

    for(auto species : mineral_species)
        ln_k[ireaction++] = lnEquilibriumConstant(species, T, P);

    return ln_k;
}

auto Phreeqx::Impl::lnActivityCoefficientsAqueousSpecies() -> Vector
{
    // Define some auxiliary variables
    const unsigned num_aqueous = aqueous_species.size();
    const double ln_10 = std::log(10.0);

    // The vector of natural log of activity coefficients of the aqueous species
    Vector ln_g(num_aqueous);

    // Update the activity coefficients of the aqueous species
    if(phreeqc.pitzer_model || phreeqc.sit_model)
    {
        // Calculate the activity coefficients using either Pitzer or SIT models
        if(phreeqc.pitzer_model)
            phreeqc.pitzer();
        else phreeqc.sit();

        // Collect the updated activity coefficients
        unsigned ispecies = 0;
        for(auto species : aqueous_species)
            ln_g[ispecies++] = species->lg_pitzer * ln_10;
    }
    else
    {
        // Calculate the activity coefficients using conventional Phreeqc models
        phreeqc.gammas(phreeqc.mu_x);

        // Collect the updated activity coefficients
        unsigned ispecies = 0;
        for(auto species : aqueous_species)
            ln_g[ispecies++] = species->lg * ln_10;
    }

    // Calculate the natural log of the activities of the aqueous species
    Vector ln_a(num_aqueous);
    for(unsigned i = 0; i < num_aqueous; ++i)
        ln_a[i] = ln_g[i] + aqueous_species[i]->lm*ln_10;

    // Calculate the total moles in the aqueous phase
    double sum = 0.0;
    for(auto species : aqueous_species)
        sum += species->moles;

    // Calculate the molar fraction of H2O
    const double nH2O = aqueous_species[iH2O]->moles;
    const double xH2O = nH2O/sum;

    // Calculate the activity of water
    if(phreeqc.pitzer_model || phreeqc.sit_model)
        ln_a[iH2O] = std::log(phreeqc.AW);
    else
        ln_a[iH2O] = std::log(xH2O);

    return ln_a;
}

auto Phreeqx::Impl::lnActivitiesAqueousSpecies() -> Vector
{
    // Define some auxiliary variables
    const unsigned num_aqueous = aqueous_species.size();
    const double ln_10 = std::log(10.0);

    // The vector of natural log of activity coefficients of the aqueous species
    Vector ln_g(num_aqueous);

    // Update the activity coefficients of the aqueous species
    if(phreeqc.pitzer_model || phreeqc.sit_model)
    {
        // Calculate the activity coefficients using either Pitzer or SIT models
        if(phreeqc.pitzer_model)
            phreeqc.pitzer();
        else phreeqc.sit();

        // Collect the updated activity coefficients
        unsigned ispecies = 0;
        for(auto species : aqueous_species)
            ln_g[ispecies++] = species->lg_pitzer * ln_10;
    }
    else
    {
        // Calculate the activity coefficients using conventional Phreeqc models
        phreeqc.gammas(phreeqc.mu_x);

        // Collect the updated activity coefficients
        unsigned ispecies = 0;
        for(auto species : aqueous_species)
            ln_g[ispecies++] = species->lg * ln_10;
    }

    // Calculate the natural log of the activities of the aqueous species
    Vector ln_a(num_aqueous);
    for(unsigned i = 0; i < num_aqueous; ++i)
        ln_a[i] = ln_g[i] + aqueous_species[i]->lm*ln_10;

    // Calculate the total moles in the aqueous phase
    double sum = 0.0;
    for(auto species : aqueous_species)
        sum += species->moles;

    // Calculate the molar fraction of H2O
    const double nH2O = aqueous_species[iH2O]->moles;
    const double xH2O = nH2O/sum;

    // Calculate the activity of water
    if(phreeqc.pitzer_model || phreeqc.sit_model)
        ln_a[iH2O] = std::log(phreeqc.AW);
    else
        ln_a[iH2O] = std::log(xH2O);

    return ln_a;
}

auto Phreeqx::Impl::lnActivitiesGaseousSpecies() -> Vector
{
    const double T = temperature();
    const double P = pressure();
    const double Patm = P * pascal_to_atm;
    const double Pbar = P * pascal_to_bar;

    phreeqc.calc_PR(gaseous_species, Patm, T, 0.0);

    const unsigned num_gaseous_species = gaseous_species.size();

    Vector ln_a(num_gaseous_species);
    for(unsigned i = 0; i < num_gaseous_species; ++i)
    {
        const double x = gaseous_species[i]->fraction_x; // the molar fraction of the gas
        const double phi = gaseous_species[i]->pr_phi;   // the fugacity coefficient of the gas
        ln_a[i] = std::log(x * phi * Pbar);
    }

    return ln_a;
}

auto Phreeqx::Impl::lnActivitiesMineralSpecies() -> Vector
{
    const unsigned num_mineral = mineral_species.size();
    return Vector::Zero(num_mineral);;
}

auto Phreeqx::Impl::lnActivities() -> Vector
{
    Vector ln_a(numSpecies());

    ln_a << lnActivitiesAqueousSpecies(),
            lnActivitiesGaseousSpecies(),
            lnActivitiesMineralSpecies();

    for(unsigned i = 0; i < numSpecies(); ++i)
        if(speciesAmount(i) == 0.0)
            ln_a[i] = 0.0;

    return ln_a;
}

auto Phreeqx::Impl::standardMolarGibbsEnergies() -> Vector
{
    /// The universal gas constant (in units of J/(mol*K))
    const double R = universalGasConstant;
    const double T = temperature();

    // Calculate the natural log of the equilibrium constants
    Vector ln_k = lnEquilibriumConstants();

    // Use the SVD decomposition of the stoichiometric matrix to calculate `u0`
    Vector u0 = svd.solve(ln_k);
    u0 *= -R*T;

    return u0;
}

auto Phreeqx::Impl::standardMolarVolumesAqueousSpecies() -> Vector
{
    const double Tc = temperature() - 273.15;
    const double Patm = pressure() * pascal_to_atm;
    const unsigned size = aqueous_species.size();
    phreeqc.calc_vm(Tc, Patm);
    Vector v(size);
    for(unsigned i = 0; i < size; ++i)
        v[i] = aqueous_species[i]->logk[vm_tc] * cm3_to_m3;
    return v;
}

auto Phreeqx::Impl::standardMolarVolumesGaseousSpecies() -> Vector
{
    const double R = universalGasConstant;
    const double T = temperature();
    const double P = pressure();
    const unsigned size = gaseous_species.size();
    return R*T/P * ones(size);
}

auto Phreeqx::Impl::standardMolarVolumesMineralSpecies() -> Vector
{
    const unsigned size = mineral_species.size();
    Vector v(size);
    for(unsigned i = 0; i < mineral_species.size(); ++i)
        v[i] = mineral_species[i]->logk[vm0] * cm3_to_m3;
    return v;
}

auto Phreeqx::Impl::standardMolarVolumes() -> Vector
{
    Vector v(numSpecies());

    v << standardMolarVolumesAqueousSpecies(),
         standardMolarVolumesGaseousSpecies(),
         standardMolarVolumesMineralSpecies();

    return v;
}

auto Phreeqx::Impl::molarVolumeAqueousPhase() -> double
{
    const Vector n_aqueous = speciesAmountsAqueousSpecies();
    const double n_total = sum(n_aqueous);
    if(n_total <= 0.0) return 0.0;
    const Vector v_aqueous = standardMolarVolumesAqueousSpecies();
    return dot(v_aqueous, n_aqueous)/n_total;
}

auto Phreeqx::Impl::molarVolumeGaseousPhase() -> double
{
    const double T = temperature();
    const double P = pressure();
    const double Patm = P * pascal_to_atm;
    const Vector n_gaseous = speciesAmountsGaseousSpecies();
    const double n_total = sum(n_gaseous);
    if(n_total <= 0.0) return 0.0;
    return phreeqc.calc_PR(gaseous_species, Patm, T, 0.0);
}

auto Phreeqx::Impl::molarVolumeMineralPhases() -> Vector
{
    return standardMolarVolumesMineralSpecies();
}

auto Phreeqx::Impl::phaseMolarVolumes() -> Vector
{
    const unsigned num_phases = numPhases();

    Vector vphases(num_phases);

    if(num_phases > 0)
    {
        unsigned offset = 0;

        if(aqueous_species.size())
            vphases[offset++] = molarVolumeAqueousPhase();

        if(gaseous_species.size())
            vphases[offset++] = molarVolumeGaseousPhase();

        rows(vphases, offset, num_phases-offset) = molarVolumeMineralPhases();
    }

    return vphases;
}

Phreeqx::Phreeqx()
: pimpl(new Impl())
{}

Phreeqx::Phreeqx(std::string database, std::string script)
: pimpl(new Impl(database, script))
{}

auto Phreeqx::set(double T, double P) -> void
{
    // Set the temperature member (in units of K)
    pimpl->phreeqc.tk_x = T;

    // Set the temperature member (in units of celsius)
    pimpl->phreeqc.tc_x = T - 298.15;

    // Set the pressure member (in units of atm)
    pimpl->phreeqc.patm_x = P * pascal_to_atm;
}

auto Phreeqx::set(double T, double P, const Vector& n) -> void
{
    set(T, P);
    pimpl->setSpeciesAmounts(n);
}

auto Phreeqx::temperature() const -> double
{
    return pimpl->temperature();
}

auto Phreeqx::pressure() const -> double
{
    return pimpl->pressure();
}

auto Phreeqx::speciesAmounts() const -> Vector
{
    return pimpl->speciesAmounts();
}

auto Phreeqx::numElements() const -> unsigned
{
    return pimpl->numElements();
}

auto Phreeqx::numSpecies() const -> unsigned
{
    return pimpl->numSpecies();
}

auto Phreeqx::numPhases() const -> unsigned
{
    return pimpl->numPhases();
}

auto Phreeqx::numSpeciesInPhase(Index iphase) const -> unsigned
{
    const unsigned num_aqueous = pimpl->aqueous_species.size();
    const unsigned num_gaseous = pimpl->gaseous_species.size();

    // Ensure `index` is not out-of-bound
    Assert(iphase < numPhases(),
        "Cannot get the number of species in phase with index "
        "`" + std::to_string(iphase) + "`.", "The given index is out of range.")

    // Return the number of aqueous species
    if(iphase == 0) return num_aqueous;

    // Return the number of gaseous species if they exist
    if(iphase == 1 && num_gaseous) return num_gaseous;

    // Return 1 as all mineral phases only have one species
    return 1;
}

auto Phreeqx::elementName(Index ielement) const -> std::string
{
    return pimpl->element_names[ielement];
}

auto Phreeqx::elementMolarMass(Index ielement) const -> double
{
    return pimpl->element_molar_masses[ielement];
}

auto Phreeqx::elementStoichiometry(Index ispecies, Index ielement) const -> double
{
    return pimpl->formula_matrix(ielement, ispecies);
}

auto Phreeqx::speciesName(Index ispecies) const -> std::string
{
    return pimpl->species_names[ispecies];
}

auto Phreeqx::phaseName(Index iphase) const -> std::string
{
    return pimpl->phase_names[iphase];
}

auto Phreeqx::phaseReferenceState(Index iphase) const -> std::string
{
    // Check if phase is the aqueous phase
    if(iphase == 0) return "IdealSolution";

    // Check if phase is the gaseous phase
    if(iphase == 1 && pimpl->gaseous_species.size()) return "IdealGas";

    // Check if phase is a mineral phase
    else return "IdealSolution";
}

auto Phreeqx::standardMolarGibbsEnergies() const -> Vector
{
    return pimpl->standardMolarGibbsEnergies();
}

auto Phreeqx::standardMolarEnthalpies() const -> Vector
{

}

auto Phreeqx::standardMolarVolumes() const -> Vector
{
    return pimpl->standardMolarVolumes();
}

auto Phreeqx::standardMolarHeatCapacitiesConstP() const -> Vector
{

}

auto Phreeqx::standardMolarHeatCapacitiesConstV() const -> Vector
{

}

auto Phreeqx::lnActivityCoefficients() const -> Vector
{
    return pimpl->lnActivityCoefficients();
}

auto Phreeqx::lnActivities() const -> Vector
{
    return pimpl->lnActivities();
}

auto Phreeqx::phaseMolarVolumes() const -> Vector
{
    return pimpl->phaseMolarVolumes();
}

auto Phreeqx::phaseResidualMolarGibbsEnergies() const -> Vector
{

}

auto Phreeqx::phaseResidualMolarEnthalpies() const -> Vector
{

}

auto Phreeqx::phaseResidualMolarHeatCapacitiesConstP() const -> Vector
{

}

auto Phreeqx::phaseResidualMolarHeatCapacitiesConstV() const -> Vector
{

}




























auto Phreeqx::standardGibbsEnergies() -> Vector
{
    return pimpl->standardMolarGibbsEnergies();
}

auto Phreeqx::standardVolumes() -> Vector
{
    return pimpl->standardMolarVolumes();
}

auto Phreeqx::activities() -> Vector
{
    return pimpl->lnActivities();
}

auto Phreeqx::chemicalPotentials() -> Vector
{
    return pimpl->chemicalPotentials();
}

auto Phreeqx::phaseMolarVolumes() -> Vector
{
    return pimpl->phaseMolarVolumes();
}

auto Phreeqx::phreeqc() -> Phreeqc&
{
    return pimpl->phreeqc;
}

auto Phreeqx::phreeqc() const -> const Phreeqc&
{
    return pimpl->phreeqc;
}

auto operator<<(std::ostream& out, const Phreeqx& phreeqx) -> std::ostream&
{
    out << "------------------------------";
    out << "Aqueous Phase";
    out << "------------------------------";
    out << std::endl;

    for(unsigned i = 0; i < phreeqx.numSpecies(); ++i)
        out << phreeqx.speciesName(i) << std::endl;

    return out;
}

} // namespace Reaktoro
