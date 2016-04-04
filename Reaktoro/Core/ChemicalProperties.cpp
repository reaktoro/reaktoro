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

#include "ChemicalProperties.hpp"

// Reaktoro includes
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Common/ReactionEquation.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
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

struct ChemicalProperties::Impl
{
    /// The chemical system
    ChemicalSystem system;

    /// The number of species in the system
    Index num_species = 0;

    /// The number of phases in the system
    Index num_phases = 0;

    /// The temperature of the system (in units of K)
    Temperature T;

    /// The pressure of the system (in units of Pa)
    Pressure P;

    /// The molar amounts of the species in the system (in units of mol).
    Vector n;

    /// The results of the evaluation of the PhaseThermoModel functions of each phase.
    std::vector<PhaseThermoModelResult> tres;

    /// The results of the evaluation of the PhaseChemicalModel functions of each phase.
    std::vector<PhaseChemicalModelResult> cres;

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a Impl instance with given ChemicalSystem
    Impl(const ChemicalSystem& system)
    : system(system)
    {
        // Initialize the number of species and phases
        num_species = system.numSpecies();
        num_phases = system.numPhases();

        // Initialize the thermodynamic and chemical properties of the phases
        tres.resize(num_phases);
        cres.resize(num_phases);
    }

    /// Update the thermodynamic properties of the chemical system.
    auto update(double T_, double P_) -> void
    {
        // Update both temperature and pressure
        T = T_;
        P = P_;

        // Update the thermodynamic properties of each phase
        for(unsigned i = 0; i < num_phases; ++i)
            tres[i] = system.phase(i).thermoModel()(T_, P_);
    }

    /// Update the chemical properties of the chemical system.
    auto update(double T_, double P_, const Vector& n_) -> void
    {
        // Set temperature, pressure and composition
        T = T_;
        P = P_;
        n = n_;

        // The offset index of the first species in each phase
        Index offset = 0;

        // Update the thermodynamic and chemical properties of each phase
        for(unsigned i = 0; i < num_phases; ++i)
        {
            // The number of species in the current phase
            const Index size = system.numSpeciesInPhase(i);

            // The vector of molar amounts of the species in the current phase
            auto np = rows(n, offset, size);

            // Calculate the phase thermodynamic and chemical properties
            tres[i] = system.phase(i).thermoModel()(T_, P_);
            cres[i] = system.phase(i).chemicalModel()(T_, P_, np);

            // Update the index of the first species in the next phase
            offset += size;
        }
    }

    /// Return the molar fractions of the species.
    auto molarFractions() const -> ChemicalVector
    {
        ChemicalVector res(num_species);
        unsigned offset = 0;
        for(unsigned i = 0; i < num_phases; ++i)
        {
            const unsigned size = system.numSpeciesInPhase(i);
            const auto np = rows(n, offset, size);
            const auto xp = Reaktoro::molarFractions(np);
            res.rows(offset, offset, size, size) = xp;
            offset += size;
        }
        return res;
    }

    /// Return the ln activity coefficients of the species.
    auto lnActivityCoefficients() const -> ChemicalVector
    {
        ChemicalVector res(num_species);
        unsigned offset = 0;
        for(unsigned i = 0; i < num_phases; ++i)
        {
            const unsigned size = system.numSpeciesInPhase(i);
            res.rows(offset, offset, size, size) = cres[i].ln_activity_coefficients;
            offset += size;
        }
        return res;
    }

    /// Return the ln activity constants of the species.
    auto lnActivityConstants() const -> ThermoVector
    {
        ThermoVector res(num_species);
        unsigned offset = 0;
        for(unsigned i = 0; i < num_phases; ++i)
        {
            const unsigned size = system.numSpeciesInPhase(i);
            res.rows(offset, size) = cres[i].ln_activity_constants;
            offset += size;
        }
        return res;
    }

    /// Return the ln activities of the species.
    auto lnActivities() const -> ChemicalVector
    {
        ChemicalVector res(num_species);
        unsigned offset = 0;
        for(unsigned i = 0; i < num_phases; ++i)
        {
            const unsigned size = system.numSpeciesInPhase(i);
            res.rows(offset, offset, size, size) = cres[i].ln_activities;
            offset += size;
        }
        return res;
    }

    /// Return the chemical potentials of the species (in units of J/mol).
    auto chemicalPotentials() const -> ChemicalVector
    {
        const auto& R = universalGasConstant;
        const auto& G = standardPartialMolarGibbsEnergies();
        const auto& lna = lnActivities();
        return G + R*T*lna;
    }

    /// Return the standard partial molar Gibbs energies of the species (in units of J/mol).
    auto standardPartialMolarGibbsEnergies() const -> ThermoVector
    {
        ThermoVector res(num_species);
        unsigned offset = 0;
        for(unsigned i = 0; i < num_phases; ++i)
        {
            const unsigned size = system.numSpeciesInPhase(i);
            res.rows(offset, size) = tres[i].standard_partial_molar_gibbs_energies;
            offset += size;
        }
        return res;
    }

    /// Return the standard partial molar enthalpies of the species (in units of J/mol).
    auto standardPartialMolarEnthalpies() const -> ThermoVector
    {
        ThermoVector res(num_species);
        unsigned offset = 0;
        for(unsigned i = 0; i < num_phases; ++i)
        {
            const unsigned size = system.numSpeciesInPhase(i);
            res.rows(offset, size) = tres[i].standard_partial_molar_enthalpies;
            offset += size;
        }
        return res;
    }

    /// Return the standard partial molar volumes of the species (in units of m3/mol).
    auto standardPartialMolarVolumes() const -> ThermoVector
    {
        ThermoVector res(num_species);
        unsigned offset = 0;
        for(unsigned i = 0; i < num_phases; ++i)
        {
            const unsigned size = system.numSpeciesInPhase(i);
            res.rows(offset, size) = tres[i].standard_partial_molar_volumes;
            offset += size;
        }
        return res;
    }

    /// Return the standard partial molar entropies of the species (in units of J/(mol*K)).
    auto standardPartialMolarEntropies() const -> ThermoVector
    {
        const auto& G = standardPartialMolarGibbsEnergies();
        const auto& H = standardPartialMolarEnthalpies();
        return (H - G)/T;
    }

    /// Return the standard partial molar internal energies of the species (in units of J/mol).
    auto standardPartialMolarInternalEnergies() const -> ThermoVector
    {
        const auto& H = standardPartialMolarEnthalpies();
        const auto& V = standardPartialMolarVolumes();
        return H - P*V;
    }

    /// Return the standard partial molar Helmholtz energies of the species (in units of J/mol).
    auto standardPartialMolarHelmholtzEnergies() const -> ThermoVector
    {
        const auto& G = standardPartialMolarGibbsEnergies();
        const auto& V = standardPartialMolarVolumes();
        return G - P*V;
    }

    /// Return the standard partial molar isobaric heat capacities of the species (in units of J/(mol*K)).
    auto standardPartialMolarHeatCapacitiesConstP() const -> ThermoVector
    {
        ThermoVector res(num_species);
        unsigned offset = 0;
        for(unsigned i = 0; i < num_phases; ++i)
        {
            const unsigned size = system.numSpeciesInPhase(i);
            res.rows(offset, size) = tres[i].standard_partial_molar_heat_capacities_cp;
            offset += size;
        }
        return res;
    }

    /// Return the standard partial molar isochoric heat capacities of the species (in units of J/(mol*K)).
    auto standardPartialMolarHeatCapacitiesConstV() const -> ThermoVector
    {
        ThermoVector res(num_species);
        unsigned offset = 0;
        for(unsigned i = 0; i < num_phases; ++i)
        {
            const unsigned size = system.numSpeciesInPhase(i);
            res.rows(offset, size) = tres[i].standard_partial_molar_heat_capacities_cv;
            offset += size;
        }
        return res;
    }

    /// Return the molar Gibbs energies of the phases (in units of J/mol).
    auto phaseMolarGibbsEnergies() const -> ChemicalVector
    {
        ChemicalVector res(num_phases, num_species);
        unsigned offset = 0;
        for(unsigned i = 0; i < num_phases; ++i)
        {
            const unsigned size = system.numSpeciesInPhase(i);
            const auto np = rows(n, offset, size);
            const auto xp = Reaktoro::molarFractions(np);
            res.row(i, offset, size) = sum(xp % tres[i].standard_partial_molar_gibbs_energies);
            res.row(i, offset, size) += cres[i].residual_molar_gibbs_energy;
            offset += size;
        }
        return res;
    }

    /// Return the molar enthalpies of the phases (in units of J/mol).
    auto phaseMolarEnthalpies() const -> ChemicalVector
    {
        ChemicalVector res(num_phases, num_species);
        unsigned offset = 0;
        for(unsigned i = 0; i < num_phases; ++i)
        {
            const unsigned size = system.numSpeciesInPhase(i);
            const auto np = rows(n, offset, size);
            const auto xp = Reaktoro::molarFractions(np);
            res.row(i, offset, size) = sum(xp % tres[i].standard_partial_molar_enthalpies);
            res.row(i, offset, size) += cres[i].residual_molar_enthalpy;
            offset += size;
        }
        return res;
    }

    /// Return the molar volumes of the phases (in units of m3/mol).
    auto phaseMolarVolumes() const -> ChemicalVector
    {
        ChemicalVector res(num_phases, num_species);
        unsigned offset = 0;
        for(unsigned i = 0; i < num_phases; ++i)
        {
            const unsigned size = system.numSpeciesInPhase(i);
            if(cres[i].molar_volume.val > 0.0)
                res.row(i, offset, size) = cres[i].molar_volume;
            else
            {
                const auto np = rows(n, offset, size);
                const auto xp = Reaktoro::molarFractions(np);
                res.row(i, offset, size) = sum(xp % tres[i].standard_partial_molar_volumes);
            }

            offset += size;
        }
        return res;
    }

    /// Return the molar entropies of the phases (in units of J/(mol*K)).
    auto phaseMolarEntropies() const -> ChemicalVector
    {
        const auto& G = phaseMolarGibbsEnergies();
        const auto& H = phaseMolarEnthalpies();
        return (H - G)/T;
    }

    /// Return the molar internal energies of the phases (in units of J/mol).
    auto phaseMolarInternalEnergies() const -> ChemicalVector
    {
        const auto& H = phaseMolarEnthalpies();
        const auto& V = phaseMolarVolumes();
        return H - P*V;
    }

    /// Return the molar Helmholtz energies of the phases (in units of J/mol).
    auto phaseMolarHelmholtzEnergies() const -> ChemicalVector
    {
        const auto& G = phaseMolarGibbsEnergies();
        const auto& V = phaseMolarVolumes();
        return G - P*V;
    }

    /// Return the molar isobaric heat capacities of the phases (in units of J/(mol*K)).
    auto phaseMolarHeatCapacitiesConstP() const -> ChemicalVector
    {
        ChemicalVector res(num_phases, num_species);
        unsigned offset = 0;
        for(unsigned i = 0; i < num_phases; ++i)
        {
            const unsigned size = system.numSpeciesInPhase(i);
            const auto np = rows(n, offset, size);
            const auto xp = Reaktoro::molarFractions(np);
            res.row(i, offset, size) = sum(xp % tres[i].standard_partial_molar_heat_capacities_cp);
            res.row(i, offset, size) += cres[i].residual_molar_heat_capacity_cp;
            offset += size;
        }
        return res;
    }

    /// Return the molar isochoric heat capacities of the phases (in units of J/(mol*K)).
    auto phaseMolarHeatCapacitiesConstV() const -> ChemicalVector
    {
        ChemicalVector res(num_phases, num_species);
        unsigned offset = 0;
        for(unsigned i = 0; i < num_phases; ++i)
        {
            const unsigned size = system.numSpeciesInPhase(i);
            const auto np = rows(n, offset, size);
            const auto xp = Reaktoro::molarFractions(np);
            res.row(i, offset, size) = sum(xp % tres[i].standard_partial_molar_heat_capacities_cv);
            res.row(i, offset, size) += cres[i].residual_molar_heat_capacity_cv;
            offset += size;
        }
        return res;
    }

    /// Return the specific Gibbs energies of the phases (in units of J/kg).
    auto phaseSpecificGibbsEnergies() const -> ChemicalVector
    {
        return phaseAmounts()/phaseMasses() % phaseMolarGibbsEnergies();
    }

    /// Return the specific enthalpies of the phases (in units of J/kg).
    auto phaseSpecificEnthalpies() const -> ChemicalVector
    {
        return phaseAmounts()/phaseMasses() % phaseMolarEnthalpies();
    }

    /// Return the specific volumes of the phases (in units of m3/kg).
    auto phaseSpecificVolumes() const -> ChemicalVector
    {
        return phaseAmounts()/phaseMasses() % phaseMolarVolumes();
    }

    /// Return the specific entropies of the phases (in units of J/(kg*K)).
    auto phaseSpecificEntropies() const -> ChemicalVector
    {
        return phaseAmounts()/phaseMasses() % phaseMolarEntropies();
    }

    /// Return the specific internal energies of the phases (in units of J/kg).
    auto phaseSpecificInternalEnergies() const -> ChemicalVector
    {
        return phaseAmounts()/phaseMasses() % phaseMolarInternalEnergies();
    }

    /// Return the specific Helmholtz energies of the phases (in units of J/kg).
    auto phaseSpecificHelmholtzEnergies() const -> ChemicalVector
    {
        return phaseAmounts()/phaseMasses() % phaseMolarHelmholtzEnergies();
    }

    /// Return the specific isobaric heat capacities of the phases (in units of J/(kg*K)).
    auto phaseSpecificHeatCapacitiesConstP() const -> ChemicalVector
    {
        return phaseAmounts()/phaseMasses() % phaseMolarHeatCapacitiesConstP();
    }

    /// Return the specific isochoric heat capacities of the phases (in units of J/(kg*K)).
    auto phaseSpecificHeatCapacitiesConstV() const -> ChemicalVector
    {
        return phaseAmounts()/phaseMasses() % phaseMolarHeatCapacitiesConstV();
    }

    /// Return the densities of the phases (in units of kg/m3).
    auto phaseDensities() const -> ChemicalVector
    {
        return phaseMasses()/(phaseAmounts() % phaseMolarVolumes());
    }

    /// Return the masses of the phases (in units of kg).
    auto phaseMasses() const -> ChemicalVector
    {
        auto nc = Reaktoro::composition(n);
        auto mm = Reaktoro::molarMasses(system.species());
        ChemicalVector res(num_phases, num_species);
        unsigned offset = 0;
        for(unsigned i = 0; i < num_phases; ++i)
        {
            const unsigned size = system.numSpeciesInPhase(i);
            auto np = nc.rows(offset, offset, size, size);
            auto mmp = rows(mm, offset, size);
            res.row(i, offset, size) = sum(mmp % np);
            offset += size;
        }
        return res;
    }

    /// Return the molar amounts of the phases (in units of mol).
    auto phaseAmounts() const -> ChemicalVector
    {
        auto nc = Reaktoro::composition(n);
        ChemicalVector res(num_phases, num_species);
        unsigned offset = 0;
        for(unsigned i = 0; i < num_phases; ++i)
        {
            const unsigned size = system.numSpeciesInPhase(i);
            auto np = nc.rows(offset, offset, size, size);
            res.row(i, offset, size) = sum(np);
            offset += size;
        }
        return res;
    }

    /// Return the volumes of the phases (in units of m3).
    auto phaseVolumes() const -> ChemicalVector
    {
        return phaseAmounts() % phaseMolarVolumes();
    }

    /// Return the volume of the system (in units of m3).
    auto volume() const -> ChemicalScalar
    {
        return sum(phaseVolumes());
    }

    /// Return the volume of a subsystem defined by some phases (in units of m3).
    auto subvolume(const Indices& iphases) const -> ChemicalScalar
    {
        return sum(phaseVolumes().rows(iphases));
    }

    /// Return the total fluid volume of the system (in units of m3).
    auto fluidVolume() const -> ChemicalScalar
    {
        const Indices iphases = system.indicesFluidPhases();
        return sum(phaseVolumes().rows(iphases));
    }

    /// Return the total solid volume of the system (in units of m3).
    auto solidVolume() const -> ChemicalScalar
    {
        const Indices iphases = system.indicesSolidPhases();
        return sum(phaseVolumes().rows(iphases));
    }
};

ChemicalProperties::ChemicalProperties()
: pimpl(new Impl())
{}

ChemicalProperties::ChemicalProperties(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

ChemicalProperties::ChemicalProperties(const ChemicalProperties& other)
: pimpl(new Impl(*other.pimpl))
{}

ChemicalProperties::~ChemicalProperties()
{}

auto ChemicalProperties::operator=(ChemicalProperties other) -> ChemicalProperties&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto ChemicalProperties::update(double T, double P) -> void
{
    pimpl->update(T, P);
}

auto ChemicalProperties::update(double T, double P, const Vector& n) -> void
{
    pimpl->update(T, P, n);
}

auto ChemicalProperties::temperature() const -> double
{
    return pimpl->T.val;
}

auto ChemicalProperties::pressure() const -> double
{
    return pimpl->P.val;
}

auto ChemicalProperties::composition() const -> const Vector&
{
    return pimpl->n;
}

auto ChemicalProperties::molarFractions() const -> ChemicalVector
{
    return pimpl->molarFractions();
}

auto ChemicalProperties::lnActivityCoefficients() const -> ChemicalVector
{
    return pimpl->lnActivityCoefficients();
}

auto ChemicalProperties::lnActivityConstants() const -> ThermoVector
{
    return pimpl->lnActivityConstants();
}

auto ChemicalProperties::lnActivities() const -> ChemicalVector
{
    return pimpl->lnActivities();
}

auto ChemicalProperties::chemicalPotentials() const -> ChemicalVector
{
    return pimpl->chemicalPotentials();
}

auto ChemicalProperties::standardPartialMolarGibbsEnergies() const -> ThermoVector
{
    return pimpl->standardPartialMolarGibbsEnergies();
}

auto ChemicalProperties::standardPartialMolarEnthalpies() const -> ThermoVector
{
    return pimpl->standardPartialMolarEnthalpies();
}

auto ChemicalProperties::standardPartialMolarVolumes() const -> ThermoVector
{
    return pimpl->standardPartialMolarVolumes();
}

auto ChemicalProperties::standardPartialMolarEntropies() const -> ThermoVector
{
    return pimpl->standardPartialMolarEntropies();
}

auto ChemicalProperties::standardPartialMolarInternalEnergies() const -> ThermoVector
{
    return pimpl->standardPartialMolarInternalEnergies();
}

auto ChemicalProperties::standardPartialMolarHelmholtzEnergies() const -> ThermoVector
{
    return pimpl->standardPartialMolarHelmholtzEnergies();
}

auto ChemicalProperties::standardPartialMolarHeatCapacitiesConstP() const -> ThermoVector
{
    return pimpl->standardPartialMolarHeatCapacitiesConstP();
}

auto ChemicalProperties::standardPartialMolarHeatCapacitiesConstV() const -> ThermoVector
{
    return pimpl->standardPartialMolarHeatCapacitiesConstV();
}

auto ChemicalProperties::phaseMolarGibbsEnergies() const -> ChemicalVector
{
    return pimpl->phaseMolarGibbsEnergies();
}

auto ChemicalProperties::phaseMolarEnthalpies() const -> ChemicalVector
{
    return pimpl->phaseMolarEnthalpies();
}

auto ChemicalProperties::phaseMolarVolumes() const -> ChemicalVector
{
    return pimpl->phaseMolarVolumes();
}

auto ChemicalProperties::phaseMolarEntropies() const -> ChemicalVector
{
    return pimpl->phaseMolarEntropies();
}

auto ChemicalProperties::phaseMolarInternalEnergies() const -> ChemicalVector
{
    return pimpl->phaseMolarInternalEnergies();
}

auto ChemicalProperties::phaseMolarHelmholtzEnergies() const -> ChemicalVector
{
    return pimpl->phaseMolarHelmholtzEnergies();
}

auto ChemicalProperties::phaseMolarHeatCapacitiesConstP() const -> ChemicalVector
{
    return pimpl->phaseMolarHeatCapacitiesConstP();
}

auto ChemicalProperties::phaseMolarHeatCapacitiesConstV() const -> ChemicalVector
{
    return pimpl->phaseMolarHeatCapacitiesConstV();
}

auto ChemicalProperties::phaseSpecificGibbsEnergies() const -> ChemicalVector
{
    return pimpl->phaseSpecificGibbsEnergies();
}

auto ChemicalProperties::phaseSpecificEnthalpies() const -> ChemicalVector
{
    return pimpl->phaseSpecificEnthalpies();
}

auto ChemicalProperties::phaseSpecificVolumes() const -> ChemicalVector
{
    return pimpl->phaseSpecificVolumes();
}

auto ChemicalProperties::phaseSpecificEntropies() const -> ChemicalVector
{
    return pimpl->phaseSpecificEntropies();
}

auto ChemicalProperties::phaseSpecificInternalEnergies() const -> ChemicalVector
{
    return pimpl->phaseSpecificInternalEnergies();
}

auto ChemicalProperties::phaseSpecificHelmholtzEnergies() const -> ChemicalVector
{
    return pimpl->phaseSpecificHelmholtzEnergies();
}

auto ChemicalProperties::phaseSpecificHeatCapacitiesConstP() const -> ChemicalVector
{
    return pimpl->phaseSpecificHeatCapacitiesConstP();
}

auto ChemicalProperties::phaseSpecificHeatCapacitiesConstV() const -> ChemicalVector
{
    return pimpl->phaseSpecificHeatCapacitiesConstV();
}

auto ChemicalProperties::phaseDensities() const -> ChemicalVector
{
    return pimpl->phaseDensities();
}

auto ChemicalProperties::phaseMasses() const -> ChemicalVector
{
    return pimpl->phaseMasses();
}

auto ChemicalProperties::phaseAmounts() const -> ChemicalVector
{
    return pimpl->phaseAmounts();
}

auto ChemicalProperties::phaseVolumes() const -> ChemicalVector
{
    return pimpl->phaseVolumes();
}

auto ChemicalProperties::volume() const -> ChemicalScalar
{
    return pimpl->volume();
}

auto ChemicalProperties::subvolume(const Indices& iphases) const -> ChemicalScalar
{
    return pimpl->subvolume(iphases);
}

auto ChemicalProperties::fluidVolume() const -> ChemicalScalar
{
    return pimpl->fluidVolume();
}

auto ChemicalProperties::solidVolume() const -> ChemicalScalar
{
    return pimpl->solidVolume();
}

auto ChemicalProperties::aqueous() const -> AqueousProperties
{
    AqueousProperties res(pimpl->system);
    res.update(pimpl->T.val, pimpl->P.val, pimpl->n);
    return res;
}

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

    /// Return the ionic strength of the system.
    auto ionicStrength() const -> ChemicalScalar
    {
        // Check there is an aqueous phase in the system
        if(!has_aqueous_phase)
            return ChemicalScalar(num_species);

        // The number of moles of water
        const double nwater = na[iwater];

        // Prepare the data for the calculation of ionic strength
        auto nc = Reaktoro::composition(na);
        auto nw = Reaktoro::amount(nwater, size, iwater);
        auto na = nc.rows(ifirst, ifirst, size, size);
        auto za = rows(charges(system.species()), ifirst, size);

        // Compute the ionic strength of the aqueous phase
        ChemicalScalar I = 0.5 * sum(na % za % za)/(nw * waterMolarMass);

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

        // The RT constant
        const ThermoScalar RT = universalGasConstant * Temperature(T);

        // The columns of the formula matrix corresponding to aqueous species
        const Matrix Aa = cols(system.formulaMatrix(), ifirst, size);

        // The ln molar amounts of aqueous species
        const Vector ln_na = log(na);

        // The weights for the weighted-LU decomposition of matrix Aa
        const Vector Wa = ln_na - min(ln_na) + 1;

        // The weighted-LU decomposition of formula matrix Aa
        LU lu(Aa, Wa);

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
