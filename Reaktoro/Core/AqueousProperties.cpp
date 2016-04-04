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

#include "AqueousProperties.hpp"

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

struct AqueousProperties::Impl
{
    /// The chemical system
    ChemicalSystem system;

    /// The boolean flag that indicates if the system has an aqueous phase
    bool has_aqueous_phase;

    /// The temperature of the system (in units of K)
    double T;

    /// The pressure of the system (in units of Pa)
    double P;

    /// The molar amounts of the aqueous species in the system (in units of mol).
    Vector na;

    /// The number of species in the system
    Index num_species = 0;

    /// The index of water species
    Index iwater;

    /// The index of hydron species
    Index ihydron;

    /// The index of charge element
    Index icharge;

    /// The index of electron species
    Index ielectron;

    /// The index of aqueous phase
    Index iaqueous;

    /// The number of aqueous species in the system
    Index size;

    /// The index of the first aqueous species in the system
    Index ifirst;

    /// The evaluation result of the aqueous PhaseThermoModel function.
    PhaseThermoModelResult tres;

    /// The evaluation result of the aqueous PhaseChemicalModel function.
    PhaseChemicalModelResult cres;

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a Impl instance with given ChemicalSystem
    Impl(const ChemicalSystem& system)
    : system(system)
    {
        // Initialize the number of species
        num_species = system.numSpecies();

        // Initialize the index of water species
        iwater = system.indexSpeciesAny(alternativeWaterNames());

        // Initialize the index of the aqueous phase
        iaqueous = system.indexPhaseWithSpecies(iwater);

        // Initialize boolean flag that indicates existence of aqueous phase
        has_aqueous_phase = iaqueous < system.numPhases();

        // Initialize other indices
        if(has_aqueous_phase)
        {
            const Phase& phase = system.phase(iaqueous);
            ihydron = phase.indexSpeciesAnyWithError(alternativeChargedSpeciesNames("H+"));
            ielectron = phase.indexSpeciesAny(alternativeChargedSpeciesNames("e-"));
            icharge = system.indexElement("Z");
            size = system.numSpeciesInPhase(iaqueous);
            ifirst = system.indexFirstSpeciesInPhase(iaqueous);
        }

        // Set the index of water to its local index in the aqueous phase
        iwater = iwater - ifirst;
    }

    /// Construct a Impl instance with given ChemicalProperties
    Impl(const ChemicalProperties& properties)
    : Impl(properties.system())
    {
        update(properties);
    }

    /// Update the aqueous properties of the chemical system.
    auto update(double T_, double P_, const Vector& n_) -> void
    {
        // Set temperature, pressure and composition of the aqueous phase
        T = T_;
        P = P_;
        na = rows(n_, ifirst, size);

        // Calculate the thermodynamic and chemical properties of aqueous phase
        tres = system.phase(iaqueous).thermoModel()(T, P);
        cres = system.phase(iaqueous).chemicalModel()(T, P, na);
    }

    /// Update the aqueous properties of the chemical system.
    auto update(const ChemicalProperties& properties) -> void
    {
        // Set temperature, pressure and composition of the aqueous phase
        T = properties.temperature();
        P = properties.pressure();
        na = rows(properties.composition(), ifirst, size);

        // Update the thermodynamic and chemical properties of aqueous phase
        tres = properties.phaseThermoModelResults()[iaqueous];
        cres = properties.phaseChemicalModelResults()[iaqueous];
    }

    /// Return the ionic strength of the system.
    auto ionicStrength() const -> ChemicalScalar
    {
        // Check there is an aqueous phase in the system
        if(!has_aqueous_phase)
            return ChemicalScalar(num_species);

        // The electrical charges of the aqueous species
        const auto za = rows(charges(system.species()), ifirst, size);

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
        rows(I.ddn, ifirst, size) = rows(I.ddn, 0, size);

        return I;
    }

    /// Return the pH of the system.
    auto pH() const -> ChemicalScalar
    {
        // Check there is an aqueous phase and hydron species in the system
        if(!has_aqueous_phase || ihydron >= size)
            return ChemicalScalar(num_species);

        // Calculate pH of the aqueous phase
        ChemicalScalar pH = -cres.ln_activities[ihydron]/ln_10;

        // Resize the derivative vector from number of aqueous species to total number of species
        pH.ddn.conservativeResize(num_species);
        rows(pH.ddn, ifirst, size) = rows(pH.ddn, 0, size);

        return pH;
    }

    /// Return the pe of the system.
    auto pE() const -> ChemicalScalar
    {
        // Check there is an aqueous phase in the system
        if(!has_aqueous_phase)
            return ChemicalScalar(num_species);

        // The ln molar amounts of aqueous species
        const Vector ln_na = log(na);

        // The columns of the formula matrix corresponding to aqueous species
        const Matrix Aa = cols(system.formulaMatrix(), ifirst, size);

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
        if(ielectron < size)
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
        rows(pe.ddn, ifirst, size) = rows(pe.ddn, 0, size);

        return pe;
    }

    /// Return the pe of the system.
    auto pE(const ReactionEquation& reaction) const -> ChemicalScalar
    {
        // Check there is an aqueous phase in the system
        if(!has_aqueous_phase)
            return ChemicalScalar(num_species);

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
        if(ielectron < size)
            G0_eminus = tres.standard_partial_molar_gibbs_energies[ielectron];

        // The pe of the system
        ChemicalScalar pe(size);

        // Loop over all species in the reaction
        for(auto pair : reaction.equation())
        {
            // Skip if current species is either e- or e[-]
            if(pair.first == "e-" || pair.first == "e[-]")
                continue;

            // Find the local index of the current species and its phase index
            const Index ispecies = system.phase(iaqueous).indexSpeciesWithError(pair.first);

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
        rows(pe.ddn, ifirst, size) = rows(pe.ddn, 0, size);

        return pe;
    }

    /// Return the reduction potential of the system (in units of V).
    auto Eh() const -> ChemicalScalar
    {
        const auto RT = universalGasConstant * Temperature(T);
        const auto F = faradayConstant;
        return ln_10*RT/F*pE();
    }

    /// Return the reduction  potential of the system calculated using a given half reaction (in units of V).
    auto Eh(std::string reaction) const -> ChemicalScalar
    {
        const auto RT = universalGasConstant * Temperature(T);
        const auto F = faradayConstant;
        return ln_10*RT/F*pE(reaction);
    }
};

AqueousProperties::AqueousProperties()
: pimpl(new Impl())
{}

AqueousProperties::AqueousProperties(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

AqueousProperties::AqueousProperties(const ChemicalProperties& properties)
: pimpl(new Impl(properties))
{}

AqueousProperties::AqueousProperties(const AqueousProperties& other)
: pimpl(new Impl(*other.pimpl))
{}

AqueousProperties::~AqueousProperties()
{}

auto AqueousProperties::operator=(AqueousProperties other) -> AqueousProperties&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto AqueousProperties::update(double T, double P, const Vector& n) -> void
{
    pimpl->update(T, P, n);
}

auto AqueousProperties::update(const ChemicalProperties& properties) -> void
{
    pimpl->update(properties);
}

auto AqueousProperties::temperature() const -> double
{
    return pimpl->T;
}

auto AqueousProperties::pressure() const -> double
{
    return pimpl->P;
}

auto AqueousProperties::composition() const -> const Vector&
{
    return pimpl->na;
}

auto AqueousProperties::ionicStrength() const -> ChemicalScalar
{
    return pimpl->ionicStrength();
}

auto AqueousProperties::pH() const -> ChemicalScalar
{
    return pimpl->pH();
}

auto AqueousProperties::pE() const -> ChemicalScalar
{
    return pimpl->pE();
}

auto AqueousProperties::pE(std::string reaction) const -> ChemicalScalar
{
    return pimpl->pE(reaction);
}

auto AqueousProperties::Eh() const -> ChemicalScalar
{
    return pimpl->Eh();
}

auto AqueousProperties::Eh(std::string reaction) const -> ChemicalScalar
{
    return pimpl->Eh(reaction);
}

} // namespace Reaktoro
