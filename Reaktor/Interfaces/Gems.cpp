// Reaktor is a C++ library for computational reaction modelling.
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

#include "Gems.hpp"

// C++ includes
#include <map>
#include <set>

// Gems includes
#define IPMGEMPLUGIN
#define NOPARTICLEARRAY
#include <Reaktor/gems3k/node.h>

// Reaktor includes
#include <Reaktor/Common/TimeUtils.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/ChemicalState.hpp>

namespace Reaktor {
namespace {

auto originalSpeciesName(const Gems& gems, unsigned index) -> std::string
{
    return gems.node().pCSD()->DCNL[index];
}

auto uniqueSpeciesNames(const Gems& gems) -> std::vector<std::string>
{
    std::map<std::string, std::set<std::string>> species_names_in_phase;

    unsigned offset = 0;
    for(unsigned iphase = 0; iphase < gems.numPhases(); ++iphase)
    {
        const std::string phase_name = gems.phaseName(iphase);
        for(unsigned i = 0; i < gems.numSpeciesInPhase(iphase); ++i)
            species_names_in_phase[phase_name].insert(originalSpeciesName(gems, offset + i));
        offset += gems.numSpeciesInPhase(iphase);
    }

    std::map<std::string, std::set<std::string>> species_found_in_phases;
    for(unsigned i = 0; i < gems.numSpecies(); ++i)
    {
        std::string species_name = originalSpeciesName(gems, i);
        for(const auto& pair : species_names_in_phase)
        {
            const auto& phase_name = pair.first;
            const auto& species_names_in_this_phase = pair.second;
            if(species_names_in_this_phase.count(species_name))
                species_found_in_phases[species_name].insert(phase_name);
        }
    }

    std::vector<std::string> species_names;
    species_names.reserve(gems.numSpecies());

    for(unsigned i = 0; i < gems.numSpecies(); ++i)
    {
        std::string species_name = originalSpeciesName(gems, i);

        if(species_found_in_phases[species_name].size() == 1)
            species_names.push_back(species_name);
        else
        {
            const Index iphase = gems.phaseIndexWithSpecies(i);
            const std::string phase_name = gems.phaseName(iphase);
            if(gems.numSpeciesInPhase(iphase) == 1)
                species_names.push_back(species_name);
            else
                species_names.push_back(species_name + "(" + phase_name + ")");
        }
    }

    return species_names;
}

} // namespace

struct Gems::Impl
{
    /// The TNode instance from Gems
    TNode node;

    /// The elapsed time of the equilibrate method (in units of s)
    double elapsed_time = 0;

    /// The options for Gems
    GemsOptions options;

    /// The unique names of the species
    std::vector<std::string> species_names;
};

Gems::Gems()
: pimpl(new Impl())
{}

Gems::Gems(std::string filename)
: pimpl(new Impl())
{
    // Initialise the GEMS `node` member
    if(pimpl->node.GEM_init(filename.c_str()))
        throw std::runtime_error("Error reading the Gems chemical system specification file.");

    // Initialise the unique names of the species
    pimpl->species_names = uniqueSpeciesNames(*this);

    //------------------------------------------------------------------------------------------------------
    // The following parameters in GEMS have to be set to extremely small values to ensure that
    // small molar amounts do not interfere with activity coefficient and chemical potential calculations
    //------------------------------------------------------------------------------------------------------
    // Reset the cutoff minimum amount of stable phase in GEMS (default: 1e-20)
    pimpl->node.pActiv()->GetActivityDataPtr()->DSM = 1e-300;

    // Set the cutoff mole amount of water-solvent for aqueous phase elimination in GEMS (default: 1e-13)
    pimpl->node.pActiv()->GetActivityDataPtr()->XwMinM = 1e-300;

    // Set the cutoff mole amount of solid sorbent for sorption phase elimination (default: 1e-13)
    pimpl->node.pActiv()->GetActivityDataPtr()->ScMinM = 1e-300;

    // Set the cutoff mole amount for elimination of DC (species) in multi-component phase (default: 1e-33)
    pimpl->node.pActiv()->GetActivityDataPtr()->DcMinM = 1e-300;

    // Set the cutoff mole amount for elimination of solution phases other than aqueous (default: 1e-20)
    pimpl->node.pActiv()->GetActivityDataPtr()->PhMinM = 1e-300;

    // Set the cutoff effective molal ionic strength for calculation of aqueous activity coefficients (default: 1e-5)
    pimpl->node.pActiv()->GetActivityDataPtr()->ICmin = 1e-300;
}

auto Gems::setTemperature(double val) -> void
{
    node().setTemperature(val);
}

auto Gems::setPressure(double val) -> void
{
    node().setPressure(val);
}

auto Gems::setSpeciesAmounts(const Vector& n) -> void
{
    node().setSpeciation(n.data());
}

auto Gems::setElementAmounts(const Vector& b) -> void
{
    // Set the molar amounts of the elements
    for(unsigned i = 0; i < numElements(); ++i)
        node().pCNode()->bIC[i] = b[i];

    // Set charge to zero
    node().pCNode()->bIC[numElements() + 1] = 0.0;
}

auto Gems::setOptions(const GemsOptions& options) -> void
{
    pimpl->options = options;
}

auto Gems::numElements() const -> unsigned
{
    return node().pCSD()->nIC - 1;
}

auto Gems::numSpecies() const -> unsigned
{
    return node().pCSD()->nDC;
}

auto Gems::numPhases() const -> unsigned
{
    return node().pCSD()->nPH;
}

auto Gems::numSpeciesInPhase(unsigned index) const -> unsigned
{
    return node().pCSD()->nDCinPH[index];
}

auto Gems::elementName(unsigned index) const -> std::string
{
    return node().pCSD()->ICNL[index];
}

auto Gems::speciesName(unsigned index) const -> std::string
{
    return pimpl->species_names[index];
}

auto Gems::phaseName(unsigned index) const -> std::string
{
    return node().pCSD()->PHNL[index];
}

auto Gems::elementIndex(std::string name) const -> unsigned
{
    unsigned index = 0;
    const unsigned size = numElements();
    for(; index < size; ++index)
        if(elementName(index) == name)
            return index;
    return size;
}

auto Gems::speciesIndex(std::string name) const -> unsigned
{
    unsigned index = 0;
    const unsigned size = numSpecies();
    for(; index < size; ++index)
        if(speciesName(index) == name)
            return index;
    return size;
}

auto Gems::phaseIndex(std::string name) const -> unsigned
{
    unsigned index = 0;
    const unsigned size = numPhases();
    for(; index < size; ++index)
        if(phaseName(index) == name)
            return index;
    return size;
}

auto Gems::phaseIndexWithSpecies(unsigned ispecies) const -> Index
{
    unsigned counter = 0;
    for(unsigned i = 0; i < numPhases(); ++i)
    {
        counter += numSpeciesInPhase(i);
        if(counter > ispecies) return i;
    }
    return numPhases();
}

auto Gems::elementAtomsInSpecies(unsigned ielement, unsigned ispecies) const -> double
{
    return node().DCaJI(ispecies, ielement);
}

auto Gems::speciesCharge(unsigned index) const -> double
{
    return elementAtomsInSpecies(numElements(), index);
}

auto Gems::elementsInSpecies(unsigned index) const -> std::map<unsigned, double>
{
    std::map<unsigned, double> elements;
    for(unsigned j = 0; j < numElements(); ++j)
        if(elementAtomsInSpecies(j, index))
            elements[j] = elementAtomsInSpecies(j, index);
    return elements;
}

auto Gems::elementMolarMass(unsigned index) const -> double
{
    return node().ICmm(index);
}

auto Gems::speciesMolarMass(unsigned index) const -> double
{
    return node().DCmm(index);
}

auto Gems::temperature() const -> double
{
    return node().Get_TK();
}

auto Gems::pressure() const -> double
{
    return node().Get_P();
}

auto Gems::elementAmounts() const -> Vector
{
    Vector b(numElements());
    for(unsigned i = 0; i < b.size(); ++i)
        b[i] = node().Get_bIC(i);
    return b;
}

auto Gems::speciesAmounts() const -> Vector
{
    Vector n(numSpecies());
    for(unsigned i = 0; i < n.size(); ++i)
        n[i] = speciesAmount(i);
    return n;
}

auto Gems::speciesAmount(unsigned index) const -> double
{
    return node().Get_nDC(index);
}

auto Gems::speciesAmountsInPhase(unsigned index) const -> Vector
{
    const unsigned size = numSpeciesInPhase(index);
    Vector np(size);
    unsigned offset = 0;
    for(unsigned i = 0; i < index; ++i)
        offset += numSpeciesInPhase(i);
    for(unsigned i = 0; i < size; ++i)
        np[i] = speciesAmount(offset + i);
    return np;
}

auto Gems::formulaMatrix() const -> Matrix
{
    const unsigned E = numElements();
    const unsigned N = numSpecies();
    Matrix A(E, N);
    for(unsigned i = 0; i < N; ++i)
        for(unsigned j = 0; j < E; ++j)
            A(j, i) = elementAtomsInSpecies(j, i);
    return A;
}

auto Gems::gibbsEnergies() -> Vector
{
    const unsigned num_species = numSpecies();
    Vector u0(num_species);
    node().updateStandardGibbsEnergies();
    ACTIVITY* ap = node().pActiv()->GetActivityDataPtr();
    for(unsigned i = 0; i < num_species; ++i)
        u0[i] = ap->tpp_G[i];
    return u0;
}

auto Gems::chemicalPotentials() -> Vector
{
    const double T = temperature();
    const double R = 8.31451; // The universal gas constant used in GEMS
    const double RT = R*T;
    const unsigned num_species = numSpecies();
    Vector u(num_species);
    node().updateStandardGibbsEnergies();
    node().initActivityCoefficients();
    node().updateConcentrations();
    node().updateActivityCoefficients();
    node().updateChemicalPotentials();
    ACTIVITY* ap = node().pActiv()->GetActivityDataPtr();
    for(unsigned i = 0; i < num_species; ++i)
        u[i] = RT*ap->F[i]; // RT factor to scale back to J/mol
    return u;
}

auto Gems::standardVolumes() -> Vector
{
    const unsigned num_species = numSpecies();
    const double cm3_to_m3 = 1e-6;
    Vector v(num_species);
    node().updateStandardVolumes();
    ACTIVITY* ap = node().pActiv()->GetActivityDataPtr();
    for(unsigned i = 0; i < num_species; ++i)
        v[i] = ap->Vol[i] * cm3_to_m3;
    return v;
}

auto Gems::equilibrate() -> void
{
    Time start = time();
    node().pCNode()->NodeStatusCH =
        pimpl->options.warmstart ? NEED_GEM_SIA : NEED_GEM_AIA;
    node().GEM_run(false);
    pimpl->elapsed_time = elapsed(start);
}

auto Gems::converged() const -> bool
{
    const auto status = node().pCNode()->NodeStatusCH;
    return status == OK_GEM_AIA or status == OK_GEM_SIA;
}

auto Gems::numIterations() const -> unsigned
{
    return node().pCNode()->IterDone;
}

auto Gems::elapsedTime() const -> double
{
    return pimpl->elapsed_time;
}

auto Gems::node() -> TNode&
{
    return pimpl->node;
}

auto Gems::node() const -> const TNode&
{
    return pimpl->node;
}

namespace helper {

auto createElement(const Gems& gems, unsigned ielement) -> Element
{
    ElementData data;
    data.name = gems.elementName(ielement);
    data.molar_mass = gems.elementMolarMass(ielement);
    return Element(data);
}

auto createSpecies(const Gems& gems, unsigned ispecies) -> Species
{
    SpeciesData data;
    data.name = gems.speciesName(ispecies);
    data.molar_mass = gems.speciesMolarMass(ispecies);
    data.charge = gems.speciesCharge(ispecies);
    data.formula = data.name;
    for(auto pair : gems.elementsInSpecies(ispecies))
    {
        data.elements.push_back(createElement(gems, pair.first));
        data.atoms.push_back(pair.second);
    }
    return Species(data);
}

auto createPhase(const Gems& gems, unsigned iphase) -> Phase
{
    PhaseData data;
    data.name = gems.phaseName(iphase);
    for(unsigned ispecies = 0; ispecies < gems.numSpecies(); ++ispecies)
        data.species.push_back(createSpecies(gems, ispecies));
    return Phase(data);
}

auto createPhases(const Gems& gems) -> std::vector<Phase>
{
    std::vector<Phase> phases;
    unsigned offset = 0;
    for(unsigned iphase = 0; iphase < gems.numPhases(); ++iphase)
    {
        PhaseData data;
        data.name = gems.phaseName(iphase);
        for(unsigned ispecies = offset; ispecies < offset + gems.numSpeciesInPhase(iphase); ++ispecies)
            data.species.push_back(createSpecies(gems, ispecies));
        phases.push_back(Phase(data));
        offset += gems.numSpeciesInPhase(iphase);
    }
    return phases;
}

} // namespace helper

Gems::operator ChemicalSystem() const
{
    Gems gems = *this;

    const unsigned num_species = gems.numSpecies();

    const Vector zero_vec = zeros(num_species);
    const Matrix zero_mat = zeros(num_species, num_species);

    ChemicalSystemData data;

    data.phases = helper::createPhases(*this);

    data.standard_gibbs_energies = [=](double T, double P) mutable -> ThermoVector
    {
        gems.setTemperature(T);
        gems.setPressure(P);
        return ThermoVector(gems.gibbsEnergies(), zero_vec, zero_vec);
    };

    data.standard_volumes = [=](double T, double P) mutable -> ThermoVector
    {
        gems.setTemperature(T);
        gems.setPressure(P);
        return ThermoVector(gems.standardVolumes(), zero_vec, zero_vec);
    };

    data.chemical_potentials = [=](double T, double P, const Vector& n) mutable -> ChemicalVector
    {
        gems.setTemperature(T);
        gems.setPressure(P);
        gems.setSpeciesAmounts(n);
        return ChemicalVector(gems.chemicalPotentials(), zero_vec, zero_vec, zero_mat);
    };

    data.ln_activities = [=](double T, double P, const Vector& n) mutable -> ChemicalVector
    {
        const double R = 8.31451;
        gems.setTemperature(T);
        gems.setPressure(P);
        gems.setSpeciesAmounts(n);
        Vector g = gems.gibbsEnergies();
        Vector u = gems.chemicalPotentials();
        Vector ln_a = (u - g)/(R*T);
        return ChemicalVector(ln_a, zero_vec, zero_vec, zero_mat);
    };

    return ChemicalSystem(data);
}

Gems::operator ChemicalState() const
{
    ChemicalSystem system = *this;
    ChemicalState state(system);
    state.setTemperature(temperature());
    state.setPressure(pressure());
    state.setSpeciesAmounts(speciesAmounts());
    return state;
}

} // namespace Reaktor
