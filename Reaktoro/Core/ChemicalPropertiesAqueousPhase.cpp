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

#include "ChemicalPropertiesAqueousPhase.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Common/ReactionEquation.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Math/LU.hpp>
#include <Reaktoro/Thermodynamics/Models/PhaseChemicalModel.hpp>
#include <Reaktoro/Thermodynamics/Models/PhaseThermoModel.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {
namespace {

// The value of ln(10)
const double ln_10 = std::log(10.0);

} // namespace

struct ChemicalPropertiesAqueousPhase::Impl
{
    /// The chemical properties instance
    ChemicalProperties properties;

    /// The chemical system
    ChemicalSystem system;

    /// The boolean flag that indicates if the system has an aqueous phase
    bool has_aqueous_phase;

    /// The number of species in the system
    Index num_species = 0;

    /// The number of aqueous species in the system
    Index num_aqueous_species = 0;

    /// The index of water species
    Index iwater;

    /// The index of hydron species
    Index ihydron;

    /// The index of charge element
    Index icharge;

    /// The index of electron species
    Index ielectron;

    /// The index of aqueous phase
    Index iaqueous_phase;

    /// The index of the first aqueous species in the system
    Index ifirst;

    /// The alkalinity reaction used to compute alkalinity.
    /// The default reaction is `Alk = Na+ + K+ + 2*Ca++ + 2*Mg++ - Cl- - 2*SO4--`.
    ReactionEquation alkalinity_reaction;

    /// The indices of the aqueous species that contribute to alkalinity
    Indices alkalinity_indices;

    /// The contribution factors of the aqueous species that contribute to alkalinity
    Vector alkalinity_factors;

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a Impl instance with given ChemicalSystem
    Impl(const ChemicalProperties& properties)
    : properties(properties), system(properties.system())
    {
        // Initialize the number of species
        num_species = system.numSpecies();

        // Initialize the index of water species
        iwater = system.indexSpeciesAny(alternativeWaterNames());

        // Initialize the index of the aqueous phase
        iaqueous_phase = system.indexPhaseWithSpecies(iwater);

        // Initialize boolean flag that indicates existence of aqueous phase
        has_aqueous_phase = iaqueous_phase < system.numPhases();

        // Initialize other indices
        if(has_aqueous_phase)
        {
            // Get a reference to the aqueous phase instance
            const Phase& phase = system.phase(iaqueous_phase);

            // Set the index of the hydron species H+ or H[+]
            ihydron = phase.indexSpeciesAny(alternativeChargedSpeciesNames("H+"));

            // Set the index of the electron species e- or e[-]
            ielectron = phase.indexSpeciesAny(alternativeChargedSpeciesNames("e-"));

            // Set the index of the element charge
            icharge = system.indexElement("Z");

            // Set the number of aqueous species
            num_aqueous_species = system.numSpeciesInPhase(iaqueous_phase);

            // Set the index of the first aqueous species
            ifirst = system.indexFirstSpeciesInPhase(iaqueous_phase);

            // Set the index of water to its local index in the aqueous phase
            iwater = iwater - ifirst;

            for(auto ion : split("Na+ K+"))
            const Index iNa = system.indexSpeciesAny(alternativeChargedSpeciesNames("Na+"));
        }
    }

    auto

    /// Return the ionic strength of the system.
    auto ionicStrength() const -> ChemicalScalar
    {
        // Check there is an aqueous phase in the system
        if(!has_aqueous_phase)
            return ChemicalScalar(num_species);

        // Get the amounts of the species
        const Vector& n = properties.composition();

        // The amounts of the aqueous species
        const Vector na = rows(n, ifirst, size);

        // The electrical charges of the aqueous species
        const Vector za = rows(charges(system.species()), ifirst, size);

        // The number of moles of water
        const double nwater = na[iwater];

        // The molar amounts of the aqueous species and its derivatives
        auto nc = Reaktoro::composition(na);

        // The molar amount of water and its derivatives
        auto nw = Reaktoro::amount(nwater, size, iwater);

        // Compute the ionic strength of the aqueous phase
        ChemicalScalar I = 0.5 * sum(nc % za % za)/(nw * waterMolarMass);

        // Resize the derivative vector from number of aqueous species to total number of species
        I.ddn.conservativeResize(num_species);
        rows(I.ddn, ifirst, num_aqueous_species) = rows(I.ddn, 0, num_aqueous_species);

        return I;
    }

    /// Return the pH of the system.
    auto pH() const -> ChemicalScalar
    {
        // Check there is an aqueous phase and hydron species in the system
        if(!has_aqueous_phase || ihydron >= num_aqueous_species)
            return ChemicalScalar(num_species);

        // Get the result of the chemical model of the aqueous phase
        const auto& cres = properties.phaseChemicalModelResults()[iaqueous_phase];

        // Calculate pH of the aqueous phase
        ChemicalScalar pH = -cres.ln_activities[ihydron]/ln_10;

        // Resize the derivative vector from number of aqueous species to total number of species
        pH.ddn.conservativeResize(num_species);
        rows(pH.ddn, ifirst, num_aqueous_species) = rows(pH.ddn, 0, num_aqueous_species);

        return pH;
    }

    /// Return the pe of the system.
    auto pE() const -> ChemicalScalar
    {
        // Check there is an aqueous phase in the system
        if(!has_aqueous_phase)
            return ChemicalScalar(num_species);

        // Get the result of the thermo and chemical models of the aqueous phase
        const auto& tres = properties.phaseThermoModelResults()[iaqueous_phase];
        const auto& cres = properties.phaseChemicalModelResults()[iaqueous_phase];

        // Get temperature of the system
        const double T = properties.temperature();

        // Get the amounts of the species
        const Vector& n = properties.composition();

        // The ln molar amounts of aqueous species
        const Vector ln_na = log(rows(n, ifirst, num_aqueous_species));

        // The columns of the formula matrix corresponding to aqueous species
        const Matrix Aa = cols(system.formulaMatrix(), ifirst, num_aqueous_species);

        // The weights for the weighted-LU decomposition of matrix Aa
        const Vector Wa = ln_na - min(ln_na) + 1;

        // The weighted-LU decomposition of formula matrix Aa
        LU lu(Aa, Wa);

        // The RT constant
        const ThermoScalar RT = universalGasConstant * Temperature(T);

        // The normalized standard chemical potentials of the aqueous species
        const ThermoVector& u0a = tres.standard_partial_molar_gibbs_energies/RT;

        // The ln activities of the aqueous species
        const ChemicalVector& ln_aa = cres.ln_activities;

        // The normalized chemical potentials of the aqueous species
        const ChemicalVector ua = u0a + ln_aa;

        // The standard chemical potential of electron species (zero if not existent in the system)
        ThermoScalar u0a_electron;
        if(ielectron < num_aqueous_species)
            u0a_electron = u0a[ielectron];

        // The dual potentials of the elements and its derivatives
        ChemicalVector y;
        y.val = lu.trsolve(ua.val);
        y.ddt = lu.trsolve(ua.ddt);
        y.ddp = lu.trsolve(ua.ddp);
        y.ddn = lu.trsolve(ua.ddn);

        // The pe of the aqueous phase
        ChemicalScalar pe(num_species);

        // The pe of the aqueous phase
        pe = (y[icharge] - u0a_electron)/ln_10;

        // Resize the derivative vector from number of aqueous species to total number of species
        pe.ddn.conservativeResize(num_species);
        rows(pe.ddn, ifirst, num_aqueous_species) = rows(pe.ddn, 0, num_aqueous_species);

        return pe;
    }

    /// Return the pe of the system.
    auto pE(const ReactionEquation& reaction) const -> ChemicalScalar
    {
        // Check there is an aqueous phase in the system
        if(!has_aqueous_phase)
            return ChemicalScalar(num_species);

        // Get the result of the thermo and chemical models of the aqueous phase
        const auto& tres = properties.phaseThermoModelResults()[iaqueous_phase];
        const auto& cres = properties.phaseChemicalModelResults()[iaqueous_phase];

        // Get the temperature of the system
        const double T = properties.temperature();

        // The RT constant
        const ThermoScalar RT = universalGasConstant * Temperature(T);

        // Find the stoichiometry of e-
        double stoichiometry_eminus = 0.0;
        if(stoichiometry_eminus == 0.0) stoichiometry_eminus = reaction.stoichiometry("e-");
        if(stoichiometry_eminus == 0.0) stoichiometry_eminus = reaction.stoichiometry("e[-]");

        // Assert the stoichiometry of e- is positive
        Assert(stoichiometry_eminus != 0.0, "Could not calculate the pe of the system.",
            "There is no `e-` or `e[-]` species in the half reaction.");

        // Find the standard chemical potential of e- (it is not zero if the
        // standard chemical potentials were obtained from log(k)'s of reactions.
        ThermoScalar G0_eminus;
        if(ielectron < num_aqueous_species)
            G0_eminus = tres.standard_partial_molar_gibbs_energies[ielectron];

        // The pe of the system
        ChemicalScalar pe(num_aqueous_species);

        // Loop over all species in the reaction
        for(auto pair : reaction.equation())
        {
            // Skip if current species is either e- or e[-]
            if(pair.first == "e-" || pair.first == "e[-]")
                continue;

            // Find the local index of the current species and its phase index
            const Index ispecies = system.phase(iaqueous_phase).indexSpeciesWithError(pair.first);

            // Get the standard chemical potential of the current species
            const ThermoScalar G0i = tres.standard_partial_molar_gibbs_energies[ispecies]/RT;

            // Get the ln activity of the current species
            const ChemicalScalar ln_ai = cres.ln_activities[ispecies];

            // Get the stoichiometry of current species
            const double stoichiometry = pair.second;

            // Update contribution
            pe -= stoichiometry * (G0i + ln_ai);
        }

        // Finalize the calculation of pe
        pe /= stoichiometry_eminus;
        pe -= G0_eminus;
        pe /= -ln_10;

        // Resize the derivative vector from number of aqueous species to total number of species
        pe.ddn.conservativeResize(num_species);
        rows(pe.ddn, ifirst, num_aqueous_species) = rows(pe.ddn, 0, num_aqueous_species);

        return pe;
    }

    /// Return the reduction potential of the system (in units of V).
    auto Eh() const -> ChemicalScalar
    {
        const double T = properties.temperature();
        const auto RT = universalGasConstant * Temperature(T);
        const auto F = faradayConstant;
        return ln_10*RT/F*pE();
    }

    /// Return the reduction  potential of the system calculated using a given half reaction (in units of V).
    auto Eh(std::string reaction) const -> ChemicalScalar
    {
        const double T = properties.temperature();
        const auto RT = universalGasConstant * Temperature(T);
        const auto F = faradayConstant;
        return ln_10*RT/F*pE(reaction);
    }
};

ChemicalPropertiesAqueousPhase::ChemicalPropertiesAqueousPhase(const ChemicalProperties& properties)
: pimpl(new Impl(properties))
{}

auto ChemicalPropertiesAqueousPhase::ionicStrength() const -> ChemicalScalar
{
    return pimpl->ionicStrength();
}

auto ChemicalPropertiesAqueousPhase::pH() const -> ChemicalScalar
{
    return pimpl->pH();
}

auto ChemicalPropertiesAqueousPhase::pE() const -> ChemicalScalar
{
    return pimpl->pE();
}

auto ChemicalPropertiesAqueousPhase::pE(std::string reaction) const -> ChemicalScalar
{
    return pimpl->pE(reaction);
}

auto ChemicalPropertiesAqueousPhase::Eh() const -> ChemicalScalar
{
    return pimpl->Eh();
}

auto ChemicalPropertiesAqueousPhase::Eh(std::string reaction) const -> ChemicalScalar
{
    return pimpl->Eh(reaction);
}

} // namespace Reaktoro
