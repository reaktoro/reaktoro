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
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Thermodynamics/Models/PhaseChemicalModel.hpp>
#include <Reaktoro/Thermodynamics/Models/PhaseThermoModel.hpp>

namespace Reaktoro {
namespace {

struct AqueousProperties
{

};

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

auto ChemicalProperties::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto ChemicalProperties::phaseThermoModelResults() const -> const std::vector<PhaseThermoModelResult>&
{
    return pimpl->tres;
}

auto ChemicalProperties::phaseChemicalModelResults() const -> const std::vector<PhaseChemicalModelResult>&
{
    return pimpl->cres;
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

} // namespace Reaktoro
