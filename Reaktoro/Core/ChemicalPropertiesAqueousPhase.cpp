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
#include <Reaktoro/Thermodynamics/Models/ChemicalModel.hpp>
#include <Reaktoro/Thermodynamics/Models/ThermoModel.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {
namespace {

struct AqueousEssentials
{
    /// The chemical properties of the system
    ChemicalProperties properties;

    /// The chemical system
    ChemicalSystem system;

    /// The boolean flag that indicates if the system has an aqueous phase
    bool has_aqueous_phase;

    /// The number of species in the system
    Index num_species = 0;

    /// The index of the aqueous phase
    Index iaqueous_phase;

    /// The number of aqueous species in the system
    Index num_aqueous_species = 0;

    /// The index of the first aqueous species in the system
    Index ifirst;

    /// The local index of the water species in the aqueous phase
    Index iwater = -1;

    /// The local index of the hydron species in the aqueous phase
    Index ihydron;

    /// The local index of the electron species in the aqueous phase
    Index ielectron;

    /// The index of the aqueous charge element
    Index icharge;

    /// The indices of the aqueous species that contribute to alkalinity
    Indices alkalinity_indices;

    /// The contribution factors of the aqueous species that contribute to alkalinity
    Vector alkalinity_factors;

    /// Construct an AqueousEssentials instance
    AqueousEssentials(const ChemicalProperties& properties)
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
            ihydron = phase.indexSpeciesAnyWithError(alternativeChargedSpeciesNames("H+"));

            // Set the index of the electron species e- or e[-]
            ielectron = phase.indexSpeciesAny(alternativeChargedSpeciesNames("e-"));

            // Set the index of the element charge
            icharge = system.indexElementWithError("Z");

            // Set the number of aqueous species
            num_aqueous_species = system.numSpeciesInPhase(iaqueous_phase);

            // Set the index of the first aqueous species
            ifirst = system.indexFirstSpeciesInPhase(iaqueous_phase);

            // Set the index of water to its local index in the aqueous phase
            iwater = iwater - ifirst;

            // Get the indices of the ions that contribute to alkalinity
            const std::map<double, std::string> ions =
            {
                {1, "Na+"}, {1, "K+"}, {2, "Ca++"}, {2, "Mg++"}, {-1, "Cl-"}, {-2, "SO4--"}
            };

            // Iterate over all ions
            for(auto pair : ions)
            {
                // Get the index of the current ion
                const Index i = phase.indexSpeciesAny(
                    alternativeChargedSpeciesNames(pair.second));

                // Store the current index if the ion is present in the aqueous phase
                if(i < num_aqueous_species)
                    alkalinity_indices.push_back(i);
            }

            // Set the alkalinity factors of the alkalinity contributors
            alkalinity_factors.resize(alkalinity_indices.size());
            auto j = 0; for(auto i : alkalinity_indices)
                alkalinity_factors[j++] = phase.species(i).charge();
        }
    }
};

// The value of ln(10)
const double ln_10 = std::log(10.0);

} // namespace

struct ChemicalPropertiesAqueousPhase::Impl : AqueousEssentials
{
    /// The temperature of the aqueous phase
    Temperature T;

    /// The pressure of the aqueous phase
    Pressure P;

    /// The amounts of the species
    Composition n;

    /// The amounts of the aqueous species
    ChemicalVectorConstRef na;

    /// The volume of the aqueous phase
    ChemicalScalar volume;

    /// The result of te aqueous thermodynamic model evaluation
    PhaseThermoModelResultConst tres;

    /// The result of te aqueous chemical model evaluation
    PhaseChemicalModelResultConst cres;

    /// Construct a Impl instance with given ChemicalSystem
    Impl(const ChemicalProperties& properties_)
    : AqueousEssentials(properties_),
      T(properties.temperature()), P(properties.pressure()),
      n(properties.composition()),
      na(rows(n, ifirst, num_aqueous_species)),
      volume(properties.phaseVolumes()[iaqueous_phase]),
      tres(properties.thermoModelResult().phaseProperties(iaqueous_phase, ifirst, num_aqueous_species)),
      cres(properties.chemicalModelResult().phaseProperties(iaqueous_phase, ifirst, num_aqueous_species))
    {}

    /// Return the ionic strength of the system.
    auto ionicStrength() const -> ChemicalScalar
    {
        // Check there is an aqueous phase in the system
        if(!has_aqueous_phase)
            return ChemicalScalar(num_species);

        // The electrical charges of the aqueous species
        const auto za = rows(charges(system.species()), ifirst, num_aqueous_species);

        // The amount of water and its derivatives
        const auto nw = na[iwater];

        // Compute the ionic strength of the aqueous phase
        return 0.5 * sum(na % za % za)/(nw * waterMolarMass);
    }

    /// Return the pH of the system.
    auto pH() const -> ChemicalScalar
    {
        // Check there is an aqueous phase and hydron species in the system
        if(!has_aqueous_phase || ihydron >= num_aqueous_species)
            return ChemicalScalar(num_species);

        // Calculate pH of the aqueous phase
        return -cres.ln_activities[ihydron]/ln_10;
    }

    /// Return the pe of the system.
    auto pE() const -> ChemicalScalar
    {
        // Check there is an aqueous phase in the system
        if(!has_aqueous_phase)
            return ChemicalScalar(num_species);

        // The ln molar amounts of aqueous species
        const Vector ln_na = log(na.val);

        // The columns of the formula matrix corresponding to aqueous species
        const auto Aa = cols(system.formulaMatrix(), ifirst, num_aqueous_species);

        // The weights for the weighted-LU decomposition of matrix Aa
        const Vector Wa = ln_na - min(ln_na) + 1;

        // The weighted-LU decomposition of formula matrix Aa
        LU lu(Aa, Wa);

        // The RT constant
        const auto RT = universalGasConstant * T;

        // The normalized standard chemical potentials of the aqueous species
        const auto u0a = tres.standard_partial_molar_gibbs_energies/RT;

        // The ln activities of the aqueous species
        const auto ln_aa = log(cres.ln_activities);

        // The normalized chemical potentials of the aqueous species
        const auto ua = u0a + ln_aa;

        // The standard chemical potential of electron species (zero if not existent in the system)
        ThermoScalar u0a_electron;
        if(ielectron < num_aqueous_species)
            u0a_electron = u0a[ielectron];

        // The dual potentials of the elements and its derivatives
        ChemicalVector y;
        y.val = lu.trsolve(ua.val);
        y.ddT = lu.trsolve(ua.ddT);
        y.ddP = lu.trsolve(ua.ddP);
        y.ddn = lu.trsolve(ua.ddn);

        // The pe of the aqueous phase
        ChemicalScalar pe(num_species);

        // The pe of the aqueous phase
        pe = (y[icharge] - u0a_electron)/ln_10;

        return pe;
    }

    /// Return the pe of the system.
    auto pE(const ReactionEquation& reaction) const -> ChemicalScalar
    {
        // Check there is an aqueous phase in the system
        if(!has_aqueous_phase)
            return ChemicalScalar(num_species);

        // The RT constant
        const auto RT = universalGasConstant * T;

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
            const auto ispecies = system.indexSpeciesWithError(pair.first);

            // Get the standard chemical potential of the current species
            const auto G0i = tres.standard_partial_molar_gibbs_energies[ispecies]/RT;

            // Get the ln activity of the current species
            const auto ln_ai = cres.ln_activities[ispecies];

            // Get the stoichiometry of current species
            const auto stoichiometry = pair.second;

            // Update contribution
            pe -= stoichiometry * (G0i + ln_ai);
        }

        // Finalize the calculation of pe
        pe /= stoichiometry_eminus;
        pe -= G0_eminus;
        pe /= -ln_10;

        return pe;
    }

    /// Return the reduction potential of the system (in units of V).
    auto Eh() const -> ChemicalScalar
    {
        const auto RT = universalGasConstant * T;
        const auto F = faradayConstant;
        return ln_10*RT/F*pE();
    }

    /// Return the reduction potential of the system calculated using a given half reaction (in units of V).
    auto Eh(std::string reaction) const -> ChemicalScalar
    {
        const auto RT = universalGasConstant * T;
        const auto F = faradayConstant;
        return ln_10*RT/F*pE(reaction);
    }

    /// Return the total alkalinity of the aqueous phase (in units of eq/L)
    auto alkalinity() const -> ChemicalScalar
    {
        const auto n_ions = rows(na, alkalinity_indices);
        const auto m3_to_liter = 1000.0;
        return sum(alkalinity_factors % n_ions)/(volume * m3_to_liter);
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

auto ChemicalPropertiesAqueousPhase::alkalinity() const -> ChemicalScalar
{
    return pimpl->alkalinity();
}

} // namespace Reaktoro
