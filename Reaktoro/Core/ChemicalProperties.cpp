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
#include <Reaktoro/Core/ChemicalPropertiesAqueousPhase.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Thermodynamics/Models/ChemicalModel.hpp>
#include <Reaktoro/Thermodynamics/Models/ThermoModel.hpp>

namespace Reaktoro {

struct ChemicalProperties::Impl
{
    /// The chemical system
    ChemicalSystem system;

    /// The number of species in the system
    Index num_species = 0;

    /// The number of phases in the system
    Index num_phases = 0;

    /// The results of the evaluation of the PhaseThermoModel functions of each phase.
    ThermoModelResult tres;

    /// The results of the evaluation of the PhaseChemicalModel functions of each phase.
    ChemicalModelResult cres;

    /// The temperature of the system (in units of K)
    Temperature T;

    /// The pressure of the system (in units of Pa)
    Pressure P;

    /// The amounts of the species in the system (in units of mol).
    Vector n;

    /// The mole fractions of the species in the system (in units of mol/mol).
    ChemicalVector x;

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a Impl instance with given ChemicalSystem
    Impl(const ChemicalSystem& system)
    : system(system),
      num_species(system.numSpecies()), num_phases(system.numPhases()),
      tres(num_species), cres(num_phases, num_species),
      T(298.15), P(1e-5), n(zeros(num_species)), x(num_species)
    {}

    /// Update the thermodynamic properties of the chemical system.
    auto update(double T_, double P_) -> void
    {
        // Update both temperature and pressure
        T = T_;
        P = P_;

        // Update the thermodynamic properties of each phase
        Index ispecies = 0;
        for(Index iphase = 0; iphase < num_phases; ++iphase)
        {
            const auto nspecies = system.numSpeciesInPhase(iphase);
            auto tp = tres.phaseProperties(ispecies, nspecies);
            system.phase(iphase).thermoModel()(tp, T, P);
            ispecies += nspecies;
        }
    }

    /// Update the chemical properties of the chemical system.
    auto update(VectorConstRef n_) -> void
    {
        // Update amounts of species
        n = n_;

        // Update the chemical properties of each phase
        Index ispecies = 0;
        for(Index iphase = 0; iphase < num_phases; ++iphase)
        {
            const auto nspecies = system.numSpeciesInPhase(iphase);
            const auto np = rows(n, ispecies, nspecies);
            const auto npc = Reaktoro::composition(np);
            auto xp = rows(x, ispecies, ispecies, nspecies, nspecies);
            auto cp = cres.phaseProperties(iphase, ispecies, nspecies);
            xp = npc/sum(npc);
            system.phase(iphase).chemicalModel()(cp, T, P, np);
            ispecies += nspecies;
        }
    }

    /// Update the chemical properties of the chemical system.
    auto update(double T_, double P_, VectorConstRef n_, const ThermoModelResult& tres_, const ChemicalModelResult& cres_) -> void
    {
        T = T_;
        P = P_;
        n = n_;
        tres = tres_;
        cres = cres_;
    }

    /// Return the molar fractions of the species.
    auto molarFractions() const -> ChemicalVector
    {
        return x;
    }

    /// Return the ln activity coefficients of the species.
    auto lnActivityCoefficients() const -> ChemicalVector
    {
        return cres.lnActivityCoefficients();
    }

    /// Return the ln activity constants of the species.
    auto lnActivityConstants() const -> ThermoVector
    {
        return tres.lnActivityConstants();
    }

    /// Return the ln activities of the species.
    auto lnActivities() const -> ChemicalVector
    {
        return cres.lnActivities();
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
        return tres.standardPartialMolarGibbsEnergies();
    }

    /// Return the standard partial molar enthalpies of the species (in units of J/mol).
    auto standardPartialMolarEnthalpies() const -> ThermoVector
    {
        return tres.standardPartialMolarEnthalpies();
    }

    /// Return the standard partial molar volumes of the species (in units of m3/mol).
    auto standardPartialMolarVolumes() const -> ThermoVector
    {
        return tres.standardPartialMolarVolumes();
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
        return tres.standardPartialMolarHeatCapacitiesConstP();
    }

    /// Return the standard partial molar isochoric heat capacities of the species (in units of J/(mol*K)).
    auto standardPartialMolarHeatCapacitiesConstV() const -> ThermoVector
    {
        return tres.standardPartialMolarHeatCapacitiesConstV();
    }

    /// Return the molar Gibbs energies of the phases (in units of J/mol).
    auto phaseMolarGibbsEnergies() const -> ChemicalVector
    {
        ChemicalVector res(num_phases, num_species);
        Index ispecies = 0;
        for(Index iphase = 0; iphase < num_phases; ++iphase)
        {
            const auto nspecies = system.numSpeciesInPhase(iphase);
            const auto np = rows(n, ispecies, nspecies);
            const auto xp = Reaktoro::molarFractions(np);
            const auto tp = tres.phaseProperties(ispecies, nspecies);
            const auto cp = cres.phaseProperties(iphase, ispecies, nspecies);
            row(res, iphase, ispecies, nspecies) = sum(xp % tp.standard_partial_molar_gibbs_energies);
            row(res, iphase, ispecies, nspecies) += cp.residual_molar_gibbs_energy;
            ispecies += nspecies;
        }
        return res;
    }

    /// Return the molar enthalpies of the phases (in units of J/mol).
    auto phaseMolarEnthalpies() const -> ChemicalVector
    {
        ChemicalVector res(num_phases, num_species);
        Index ispecies = 0;
        for(Index iphase = 0; iphase < num_phases; ++iphase)
        {
            const auto nspecies = system.numSpeciesInPhase(iphase);
            const auto np = rows(n, ispecies, nspecies);
            const auto xp = Reaktoro::molarFractions(np);
            const auto tp = tres.phaseProperties(ispecies, nspecies);
            const auto cp = cres.phaseProperties(iphase, ispecies, nspecies);
            row(res, iphase, ispecies, nspecies) = sum(xp % tp.standard_partial_molar_enthalpies);
            row(res, iphase, ispecies, nspecies) += cp.residual_molar_enthalpy;
            ispecies += nspecies;
        }
        return res;
    }

    /// Return the molar volumes of the phases (in units of m3/mol).
    auto phaseMolarVolumes() const -> ChemicalVector
    {
        ChemicalVector res(num_phases, num_species);
        Index ispecies = 0;
        for(Index iphase = 0; iphase < num_phases; ++iphase)
        {
            const auto nspecies = system.numSpeciesInPhase(iphase);
            const auto tp = tres.phaseProperties(ispecies, nspecies);
            const auto cp = cres.phaseProperties(iphase, ispecies, nspecies);
            if(cp.molar_volume > 0.0)
                row(res, iphase, ispecies, nspecies) = cp.molar_volume;
            else
            {
                const auto np = rows(n, ispecies, nspecies);
                const auto xp = Reaktoro::molarFractions(np);
                row(res, iphase, ispecies, nspecies) = sum(xp % tp.standard_partial_molar_volumes);
            }

            ispecies += nspecies;
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
        Index ispecies = 0;
        for(Index iphase = 0; iphase < num_phases; ++iphase)
        {
            const auto nspecies = system.numSpeciesInPhase(iphase);
            const auto np = rows(n, ispecies, nspecies);
            const auto xp = Reaktoro::molarFractions(np);
            const auto tp = tres.phaseProperties(ispecies, nspecies);
            const auto cp = cres.phaseProperties(iphase, ispecies, nspecies);
            row(res, iphase, ispecies, nspecies) = sum(xp % tp.standard_partial_molar_heat_capacities_cp);
            row(res, iphase, ispecies, nspecies) += cp.residual_molar_heat_capacity_cp;
            ispecies += nspecies;
        }
        return res;
    }

    /// Return the molar isochoric heat capacities of the phases (in units of J/(mol*K)).
    auto phaseMolarHeatCapacitiesConstV() const -> ChemicalVector
    {
        ChemicalVector res(num_phases, num_species);
        Index ispecies = 0;
        for(Index iphase = 0; iphase < num_phases; ++iphase)
        {
            const auto nspecies = system.numSpeciesInPhase(iphase);
            const auto np = rows(n, ispecies, nspecies);
            const auto xp = Reaktoro::molarFractions(np);
            const auto tp = tres.phaseProperties(ispecies, nspecies);
            const auto cp = cres.phaseProperties(iphase, ispecies, nspecies);
            row(res, iphase, ispecies, nspecies) = sum(xp % tp.standard_partial_molar_heat_capacities_cv);
            row(res, iphase, ispecies, nspecies) += cp.residual_molar_heat_capacity_cv;
            ispecies += nspecies;
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
        const auto nc = Reaktoro::composition(n);
        const auto mm = Reaktoro::molarMasses(system.species());
        ChemicalVector res(num_phases, num_species);
        Index ispecies = 0;
        for(Index iphase = 0; iphase < num_phases; ++iphase)
        {
            const auto nspecies = system.numSpeciesInPhase(iphase);
            const auto np = rows(nc, ispecies, ispecies, nspecies, nspecies);
            auto mmp = rows(mm, ispecies, nspecies);
            row(res, iphase, ispecies, nspecies) = sum(mmp % np);
            ispecies += nspecies;
        }
        return res;
    }

    /// Return the molar amounts of the phases (in units of mol).
    auto phaseAmounts() const -> ChemicalVector
    {
        const auto nc = Reaktoro::composition(n);
        ChemicalVector res(num_phases, num_species);
        Index ispecies = 0;
        for(Index iphase = 0; iphase < num_phases; ++iphase)
        {
            const auto nspecies = system.numSpeciesInPhase(iphase);
            const auto np = rows(nc, ispecies, ispecies, nspecies, nspecies);
            row(res, iphase, ispecies, nspecies) = sum(np);
            ispecies += nspecies;
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
        return sum(rows(phaseVolumes(), iphases));
    }

    /// Return the total fluid volume of the system (in units of m3).
    auto fluidVolume() const -> ChemicalScalar
    {
        const Indices iphases = system.indicesFluidPhases();
        return sum(rows(phaseVolumes(), iphases));
    }

    /// Return the total solid volume of the system (in units of m3).
    auto solidVolume() const -> ChemicalScalar
    {
        const Indices iphases = system.indicesSolidPhases();
        return sum(rows(phaseVolumes(), iphases));
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

auto ChemicalProperties::update(VectorConstRef n) -> void
{
    pimpl->update(n);
}

auto ChemicalProperties::update(double T, double P, VectorConstRef n) -> void
{
    pimpl->update(T, P);
    pimpl->update(n);
}

auto ChemicalProperties::update(double T, double P, VectorConstRef n, const ThermoModelResult& tres, const ChemicalModelResult& cres) -> void
{
    pimpl->update(T, P, n, tres, cres);
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

auto ChemicalProperties::thermoModelResult() const -> const ThermoModelResult&
{
    return pimpl->tres;
}

auto ChemicalProperties::chemicalModelResult() const -> const ChemicalModelResult&
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

auto ChemicalProperties::aqueous() const -> ChemicalPropertiesAqueousPhase
{
    ChemicalPropertiesAqueousPhase aqueous(*this);
    return aqueous;
}

} // namespace Reaktoro
