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

#include "Gems.hpp"

// C++ includes
#include <map>
#include <set>

// Gems includes
#define IPMGEMPLUGIN
#define NOPARTICLEARRAY
#include <gems/node.h>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/TimeUtils.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {
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
            const Index iphase = gems.indexPhaseWithSpecies(i);
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

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a default Impl instance
    Impl(std::string filename)
    {
        // Initialize the GEMS `node` member
        if(node.GEM_init(filename.c_str()))
            throw std::runtime_error("Error reading the Gems chemical system specification file.");

        //------------------------------------------------------------------------------------------------------
        // The following parameters in GEMS have to be set to extremely small values to ensure that
        // small molar amounts do not interfere with activity coefficient and chemical potential calculations
        //------------------------------------------------------------------------------------------------------
        // Reset the cutoff minimum amount of stable phase in GEMS (default: 1e-20)
        node.pActiv()->GetActivityDataPtr()->DSM = 1e-300;

        // Set the cutoff mole amount of water-solvent for aqueous phase elimination in GEMS (default: 1e-13)
        node.pActiv()->GetActivityDataPtr()->XwMinM = 1e-300;

        // Set the cutoff mole amount of solid sorbent for sorption phase elimination (default: 1e-13)
        node.pActiv()->GetActivityDataPtr()->ScMinM = 1e-300;

        // Set the cutoff mole amount for elimination of DC (species) in multi-component phase (default: 1e-33)
        node.pActiv()->GetActivityDataPtr()->DcMinM = 1e-300;

        // Set the cutoff mole amount for elimination of solution phases other than aqueous (default: 1e-20)
        node.pActiv()->GetActivityDataPtr()->PhMinM = 1e-300;

        // Set the cutoff effective molal ionic strength for calculation of aqueous activity coefficients (default: 1e-5)
        node.pActiv()->GetActivityDataPtr()->ICmin = 1e-300;
    }
};

Gems::Gems()
: pimpl(new Impl())
{}

Gems::Gems(std::string filename)
: pimpl(new Impl(filename))
{
    // Initialize the unique names of the species
    pimpl->species_names = uniqueSpeciesNames(*this);
}

Gems::~Gems()
{}

auto Gems::set(double T, double P) -> void
{
    node().setTemperature(T);
    node().setPressure(P);
}

auto Gems::set(double T, double P, const Vector& n) -> void
{
    node().setTemperature(T);
    node().setPressure(P);
    node().setSpeciation(n.data());

    node().updateStandardGibbsEnergies();
    node().initActivityCoefficients();
    node().updateConcentrations();
    node().updateActivityCoefficients();
    node().updateChemicalPotentials();
    node().updateActivities();
}

auto Gems::temperature() const -> double
{
    return node().Get_TK();
}

auto Gems::pressure() const -> double
{
    return node().Get_P();
}

auto Gems::speciesAmounts() const -> Vector
{
    Vector n(numSpecies());
    for(unsigned i = 0; i < n.size(); ++i)
        n[i] = node().Get_nDC(i);
    return n;
}

auto Gems::numElements() const -> unsigned
{
    return node().pCSD()->nIC;
}

auto Gems::numSpecies() const -> unsigned
{
    return node().pCSD()->nDC;
}

auto Gems::numPhases() const -> unsigned
{
    return node().pCSD()->nPH;
}

auto Gems::numSpeciesInPhase(Index iphase) const -> unsigned
{
    return node().pCSD()->nDCinPH[iphase];
}

auto Gems::elementName(Index ielement) const -> std::string
{
    return node().pCSD()->ICNL[ielement];
}

auto Gems::elementMolarMass(Index ielement) const -> double
{
    return node().ICmm(ielement);
}

auto Gems::elementStoichiometry(Index ispecies, Index ielement) const -> double
{
    return node().DCaJI(ispecies, ielement);
}

auto Gems::speciesName(Index ispecies) const -> std::string
{
    return pimpl->species_names[ispecies];
}

auto Gems::phaseName(Index iphase) const -> std::string
{
    return node().pCSD()->PHNL[iphase];
}

auto Gems::phaseReferenceState(Index iphase) const -> std::string
{
    // Check if the phase is a mixture of gases (PH_GASMIX) or fluid phase (PH_FLUID)
    ACTIVITY* ap = node().pActiv()->GetActivityDataPtr();
    if(ap->PHC[iphase] == 'g' || ap->PHC[iphase] == 'f') return "IdealGas";
    else return "IdealSolution";
}

auto Gems::standardMolarGibbsEnergies() const -> Vector
{
    const double T = temperature();
    const double P = pressure();
    TNode& nod = const_cast<TNode&>(node());
    Vector res(numSpecies());
    for(unsigned i = 0; i < res.size(); ++i)
        res[i] = nod.DC_G0(i, P, T, false);
    return res;
}

auto Gems::standardMolarEnthalpies() const -> Vector
{
    const double T = temperature();
    const double P = pressure();
    TNode& nod = const_cast<TNode&>(node());
    Vector res(numSpecies());
    for(unsigned i = 0; i < res.size(); ++i)
        res[i] = nod.DC_H0(i, P, T);
    return res;
}

auto Gems::standardMolarVolumes() const -> Vector
{
    const double cm3_to_m3 = 1e-6;
    const double T = temperature();
    const double P = pressure();
    TNode& nod = const_cast<TNode&>(node());
    Vector res(numSpecies());
    for(unsigned i = 0; i < res.size(); ++i)
        res[i] = nod.DC_V0(i, P, T) * cm3_to_m3;
    return res;
}

auto Gems::standardMolarHeatCapacitiesConstP() const -> Vector
{
    const double T = temperature();
    const double P = pressure();
    TNode& nod = const_cast<TNode&>(node());
    Vector res(numSpecies());
    for(unsigned i = 0; i < res.size(); ++i)
        res[i] = nod.DC_Cp0(i, P, T);
    return res;
}

auto Gems::standardMolarHeatCapacitiesConstV() const -> Vector
{
    return standardMolarHeatCapacitiesConstP();
}

auto Gems::lnActivityCoefficients() const -> Vector
{
    ACTIVITY* ap = node().pActiv()->GetActivityDataPtr();
    const unsigned nspecies = numSpecies();
    Vector res(nspecies);
    for(unsigned i = 0; i < nspecies; ++i)
        res[i] = ap->lnGam[i];
    return res;
}

auto Gems::lnActivities() const -> Vector
{
    ACTIVITY* ap = node().pActiv()->GetActivityDataPtr();
    const unsigned nspecies = numSpecies();
    Vector res(nspecies);
    for(unsigned i = 0; i < nspecies; ++i)
        res[i] = ap->lnAct[i];
    return res;
}

auto Gems::setOptions(const GemsOptions& options) -> void
{
    pimpl->options = options;
}

auto Gems::equilibrate(double T, double P, const Vector& b) -> void
{
    // Start timing
    Time start = time();

    // Set temperature and pressure
    set(T, P);

    // Set the molar amounts of the elements
    for(unsigned i = 0; i < numElements(); ++i)
        node().pCNode()->bIC[i] = b[i];

    // Solve the equilibrium problem with gems
    node().pCNode()->NodeStatusCH =
        pimpl->options.warmstart ? NEED_GEM_SIA : NEED_GEM_AIA;
    node().GEM_run(false);

    // Finish timing
    pimpl->elapsed_time = elapsed(start);
}

auto Gems::converged() const -> bool
{
    const auto status = node().pCNode()->NodeStatusCH;
    return status == OK_GEM_AIA || status == OK_GEM_SIA;
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

} // namespace Reaktoro
