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

#include "ChemicalProperty.hpp"

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
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {
namespace {

auto indexWaterSpecies(const ChemicalSystem& system) -> Index
{
    return system.indexSpeciesAny(alternativeWaterNames());
}

auto indexAqueousPhase(const ChemicalSystem& system) -> Index
{
    return system.indexPhaseWithSpecies(indexWaterSpecies(system));
}

const double ln_10 = std::log(10.0);

} // namespace

auto ChemicalProperty::ionicStrength(const ChemicalSystem& system) -> ChemicalPropertyFunction
{
    // The index of the aqueous phase
    const Index iaqueousphase = indexAqueousPhase(system);

    // The number of species in the system
    const Index num_species = system.numSpecies();

    // Check if there is an aqueous phase in the system
    if(iaqueousphase >= system.numPhases())
        return [=](const ChemicalProperties&) { return ChemicalScalar(num_species); };

    // The index of the first aqueous species
    const Index ifirst = system.indexFirstSpeciesInPhase(iaqueousphase);

    // The number of aqueous species
    const Index num_aqueous = system.numSpeciesInPhase(iaqueousphase);

    // The index of water species
    const Index iwater = indexWaterSpecies(system);

    // The electrical charges of the aqueous species
    const Vector za = rows(charges(system.species()), ifirst, num_aqueous);

    ChemicalPropertyFunction f = [=](const ChemicalProperties& properties)
    {
        const auto n = properties.composition();
        const auto na = rows(n, ifirst, num_aqueous);
        const auto nw = n[iwater];
        ChemicalScalar res = 0.5 * sum(na % za % za)/(nw * waterMolarMass);
        return res;
    };

    return f;
}

auto ChemicalProperty::pH(const ChemicalSystem& system) -> ChemicalPropertyFunction
{
    const Index iaqueousphase = indexAqueousPhase(system);
    const Index num_species = system.numSpecies();
    const Index ihydron = system.indexSpeciesAny(alternativeChargedSpeciesNames("H+"));

    // Check if there is an aqueous phase in the system
    if(iaqueousphase >= system.numPhases())
        return [=](const ChemicalProperties&) { return ChemicalScalar(num_species); };

    ChemicalPropertyFunction f = [=](const ChemicalProperties& properties)
    {
        ChemicalScalar res = -properties.lnActivities()[ihydron]/ln_10;
        return res;
    };

    return f;
}

auto ChemicalProperty::pE(const ChemicalSystem& system) -> ChemicalPropertyFunction
{
    const Index iaqueousphase = indexAqueousPhase(system);
    const Index num_aqueous = system.numSpeciesInPhase(iaqueousphase);
    const Index ifirst = system.indexFirstSpeciesInPhase(iaqueousphase);
    const Index ielectron = system.phase(iaqueousphase).indexSpeciesAny(alternativeChargedSpeciesNames("e-"));
    const Index num_species = system.numSpecies();
    const Index icharge = system.indexElementWithError("Z");

    // Check if there is an aqueous phase in the system
    if(iaqueousphase >= system.numPhases())
        return [=](const ChemicalProperties&) { return ChemicalScalar(num_species); };

    ChemicalPropertyFunction f = [=](const ChemicalProperties& properties)
    {
        // The amounts of the aqueous species
        const auto na = rows(properties.composition(), ifirst, num_aqueous);

        // The ln amounts of aqueous species
        const Vector ln_na = log(na.val);

        // The columns of the formula matrix corresponding to aqueous species
        const auto Aa = cols(system.formulaMatrix(), ifirst, num_aqueous);

        // The weights for the weighted-LU decomposition of matrix Aa
        const Vector Wa = ln_na - min(ln_na) + 1;

        // The weighted-LU decomposition of formula matrix Aa
        LU lu(Aa, Wa);

        // The RT constant
        const auto T = properties.temperature();
        const auto RT = universalGasConstant * T;

        // The normalized standard chemical potentials of the aqueous species
        const ThermoVector u0a = rows(properties.standardPartialMolarGibbsEnergies(), ifirst, num_aqueous)/RT;

        // The ln activities of the aqueous species
        const ChemicalVector ln_aa = rows(properties.lnActivities(), ifirst, num_aqueous);

        // The normalized chemical potentials of the aqueous species
        const ChemicalVector ua = u0a + ln_aa;

        // The standard chemical potential of electron species (zero if not existent in the system)
        ThermoScalar u0a_electron;
        if(ielectron < num_aqueous)
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
    };

    return f;
}

auto ChemicalProperty::pE(const ChemicalSystem& system, const ReactionEquation& reaction) -> ChemicalPropertyFunction
{
    const Index iaqueousphase = indexAqueousPhase(system);
    const Index ielectron = system.indexSpeciesAny(alternativeChargedSpeciesNames("e-"));
    const Index num_species = system.numSpecies();

    // Check if there is an aqueous phase in the system
    if(iaqueousphase >= system.numPhases())
        return [=](const ChemicalProperties&) { return ChemicalScalar(num_species); };

    // Find the stoichiometry of e-
    double stoichiometry_eminus = 0.0;
    if(stoichiometry_eminus == 0.0) stoichiometry_eminus = reaction.stoichiometry("e-");
    if(stoichiometry_eminus == 0.0) stoichiometry_eminus = reaction.stoichiometry("e[-]");

    // Assert the stoichiometry of e- is positive
    Assert(stoichiometry_eminus != 0.0, "Could not calculate the pe of the system.",
        "There is no `e-` or `e[-]` species in the half reaction.");

    ChemicalPropertyFunction f = [=](const ChemicalProperties& properties)
    {
        // The pE of the aqueous phase
        ChemicalScalar pe(num_species);

        const auto T = properties.temperature();
        const auto RT = universalGasConstant * T;

        // Find the standard chemical potential of e- (it is not zero if the
        // standard chemical potentials were obtained from log(k)'s of reactions.
        ThermoScalar G0_eminus;
        if(ielectron < num_species)
            G0_eminus = properties.standardPartialMolarGibbsEnergies()[ielectron];

        // Loop over all species in the reaction
        for(auto pair : reaction.equation())
        {
            // Skip if current species is either e- or e[-]
            if(pair.first == "e-" || pair.first == "e[-]")
                continue;

            // Find the local index of the current species and its phase index
            const auto ispecies = system.indexSpeciesWithError(pair.first);

            // Get the standard chemical potential of the current species
            const auto G0i = properties.standardPartialMolarGibbsEnergies()[ispecies]/RT;

            // Get the ln activity of the current species
            const auto ln_ai = properties.lnActivities()[ispecies];

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
    };

    return f;
}

auto ChemicalProperty::Eh(const ChemicalSystem& system) -> ChemicalPropertyFunction
{
    auto pE = ChemicalProperty::pE(system);

    ChemicalPropertyFunction f = [=](const ChemicalProperties& properties)
    {
        const auto T = properties.temperature();
        const auto RT = universalGasConstant * T;
        const auto F = faradayConstant;
        return ln_10*RT/F*pE(properties);
    };

    return f;
}

auto ChemicalProperty::Eh(const ChemicalSystem& system, const ReactionEquation& reaction) -> ChemicalPropertyFunction
{
    auto pE = ChemicalProperty::pE(system, reaction);

    ChemicalPropertyFunction f = [=](const ChemicalProperties& properties)
    {
        const auto T = properties.temperature();
        const auto RT = universalGasConstant * T;
        const auto F = faradayConstant;
        ChemicalScalar res = ln_10*RT/F*pE(properties);
        return res;
    };

    return f;
}

auto ChemicalProperty::alkalinity(const ChemicalSystem& system) -> ChemicalPropertyFunction
{
    // The index of the aqueous phase
    const Index iaqueousphase = indexAqueousPhase(system);

    // The number of species in the system
    const Index num_species = system.numSpecies();

    // Check if there is an aqueous phase in the system
    if(iaqueousphase >= system.numPhases())
        return [=](const ChemicalProperties&) { return ChemicalScalar(num_species); };

    // The ions that contribute to alkalinity
    const std::map<double, std::string> ions = { {1, "Na+"}, {1, "K+"}, {2, "Ca++"}, {2, "Mg++"}, {-1, "Cl-"}, {-2, "SO4--"} };

    /// The indices of the aqueous species that contribute to alkalinity
    Indices alkalinity_indices;

    // Iterate over all above ions
    for(auto pair : ions)
    {
        // Get the index of the current ion
        const Index i = system.indexSpeciesAny(
            alternativeChargedSpeciesNames(pair.second));

        // Store the current index if the ion is present in the aqueous phase
        if(i < num_species)
            alkalinity_indices.push_back(i);
    }

    /// The contribution factors of the aqueous species that contribute to alkalinity
    Vector alkalinity_factors(alkalinity_indices.size());
    auto j = 0; for(auto i : alkalinity_indices)
        alkalinity_factors[j++] = system.species(i).charge();

    ChemicalScalar volume(num_species);
    ChemicalVector n;

    ChemicalPropertyFunction f = [=](const ChemicalProperties& properties) mutable
    {
        n = properties.composition();
        const auto n_ions = rows(n, alkalinity_indices);
        const auto m3_to_liter = 1000.0;
        volume = properties.phaseVolumes()[iaqueousphase];
        ChemicalScalar res = sum(alkalinity_factors % n_ions)/(volume * m3_to_liter);
        return res;
    };

    return f;
}

} // namespace Reaktoro
