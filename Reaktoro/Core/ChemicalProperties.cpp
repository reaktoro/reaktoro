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

#include "ChemicalProperties.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {

ChemicalProperties::ChemicalProperties()
{}

ChemicalProperties::ChemicalProperties(const ChemicalSystem& system)
: system(system), num_species(system.numSpecies()), num_phases(system.numPhases()),
  T(NAN), P(NAN), n(zeros(num_species)), x(num_species),
  tres(num_phases, num_species), cres(num_phases, num_species)
{}

auto ChemicalProperties::update(double T_, double P_) -> void
{
    // Update both temperature and pressure
    if(T != T_ || P != P_)
    {
        T = T_;
        P = P_;
        system.thermoModel()(tres, T, P);
    }
}

auto ChemicalProperties::update(VectorXrConstRef n_) -> void
{
    Assert(!std::isnan(T.val) && !std::isnan(P.val),
           "Cannot proceed with method ChemicalProperties::update.",
           "The temperature or pressure values are invalid (NAN). "
           "Update these properties before calling this method!")

    n = n_;
    system.chemicalModel()(cres, T, P, n);

    // Update mole fractions
    Index offset = 0;
    for(Index iphase = 0; iphase < num_phases; ++iphase)
    {
        const auto size = system.numSpeciesInPhase(iphase);
        const auto np = rows(n, offset, size);
        auto xp = rows(x, offset, offset, size, size);
        if(size == 1) {
            xp = 1.0;
        }
        else {
            const auto snp = sum(np);
            if(snp != 0.0)
                xp = np/snp;
            else
                xp = 0.0;
        }
        offset += size;
    }
}

auto ChemicalProperties::update(double T, double P, VectorXrConstRef n) -> void
{
    update(T, P);
    update(n);
}

auto ChemicalProperties::update(double T_, double P_, VectorXrConstRef n_, const ThermoModelResult& tres_, const ChemicalModelResult& cres_) -> void
{
    T = T_;
    P = P_;
    n = n_;
    tres = tres_;
    cres = cres_;
}

auto ChemicalProperties::temperature() const -> real
{
    return T;
}

auto ChemicalProperties::pressure() const -> real
{
    return P;
}

auto ChemicalProperties::composition() const -> VectorXrConstRef
{
    return n;
}

auto ChemicalProperties::thermoModelResult() const -> const ThermoModelResult&
{
    return tres;
}

auto ChemicalProperties::chemicalModelResult() const -> const ChemicalModelResult&
{
    return cres;
}

auto ChemicalProperties::moleFractions() const -> VectorXd
{
    return x;
}

auto ChemicalProperties::lnActivityCoefficients() const -> VectorXdConstRef
{
    return cres.lnActivityCoefficients();
}

auto ChemicalProperties::lnActivityConstants() const -> VectorXrConstRef
{
    return tres.lnActivityConstants();
}

auto ChemicalProperties::lnActivities() const -> VectorXdConstRef
{
    return cres.lnActivities();
}

auto ChemicalProperties::partialMolarVolumes() const -> VectorXdConstRef
{
    return cres.partialMolarVolumes();
}

auto ChemicalProperties::chemicalPotentials() const -> VectorXd
{
    const auto& R = universalGasConstant;
    const auto& G = standardPartialMolarGibbsEnergies();
    const auto& lna = lnActivities();
    return G + R*T*lna;
}

auto ChemicalProperties::standardPartialMolarGibbsEnergies() const -> VectorXrConstRef
{
    return tres.standardPartialMolarGibbsEnergies();
}

auto ChemicalProperties::standardPartialMolarEnthalpies() const -> VectorXrConstRef
{
    return tres.standardPartialMolarEnthalpies();
}

auto ChemicalProperties::standardPartialMolarVolumes() const -> VectorXrConstRef
{
    return tres.standardPartialMolarVolumes();
}

auto ChemicalProperties::standardPartialMolarEntropies() const -> VectorXr
{
    const auto& G = standardPartialMolarGibbsEnergies();
    const auto& H = standardPartialMolarEnthalpies();
    return (H - G)/T;
}

auto ChemicalProperties::standardPartialMolarInternalEnergies() const -> VectorXr
{
    const auto& H = standardPartialMolarEnthalpies();
    const auto& V = standardPartialMolarVolumes();
    return H - P*V;
}

auto ChemicalProperties::standardPartialMolarHelmholtzEnergies() const -> VectorXr
{
    const auto& G = standardPartialMolarGibbsEnergies();
    const auto& V = standardPartialMolarVolumes();
    return G - P*V;
}

auto ChemicalProperties::standardPartialMolarHeatCapacitiesConstP() const -> VectorXrConstRef
{
    return tres.standardPartialMolarHeatCapacitiesConstP();
}

auto ChemicalProperties::standardPartialMolarHeatCapacitiesConstV() const -> VectorXrConstRef
{
    return tres.standardPartialMolarHeatCapacitiesConstV();
}

auto ChemicalProperties::phaseMolarGibbsEnergies() const -> VectorXd
{
    VectorXd res(num_phases, num_species);
    Index ispecies = 0;
    for(Index iphase = 0; iphase < num_phases; ++iphase)
    {
        const auto nspecies = system.numSpeciesInPhase(iphase);
        const auto xp = rows(x, ispecies, ispecies, nspecies, nspecies);
        const auto tp = tres.phaseProperties(iphase, ispecies, nspecies);
        const auto cp = cres.phaseProperties(iphase, ispecies, nspecies);
        row(res, iphase, ispecies, nspecies) = sum(xp % tp.standard_partial_molar_gibbs_energies);
        row(res, iphase, ispecies, nspecies) += cp.residual_molar_gibbs_energy;
        ispecies += nspecies;
    }
    return res;
}

auto ChemicalProperties::phaseMolarEnthalpies() const -> VectorXd
{
    VectorXd res(num_phases, num_species);
    Index ispecies = 0;
    for(Index iphase = 0; iphase < num_phases; ++iphase)
    {
        const auto nspecies = system.numSpeciesInPhase(iphase);
        const auto xp = rows(x, ispecies, ispecies, nspecies, nspecies);
        const auto tp = tres.phaseProperties(iphase, ispecies, nspecies);
        const auto cp = cres.phaseProperties(iphase, ispecies, nspecies);
        row(res, iphase, ispecies, nspecies) = sum(xp % tp.standard_partial_molar_enthalpies);
        row(res, iphase, ispecies, nspecies) += cp.residual_molar_enthalpy;
        ispecies += nspecies;
    }
    return res;
}

auto ChemicalProperties::phaseMolarVolumes() const -> VectorXd
{
    VectorXd res(num_phases, num_species);
    Index ispecies = 0;
    for(Index iphase = 0; iphase < num_phases; ++iphase)
    {
        const auto nspecies = system.numSpeciesInPhase(iphase);
        const auto tp = tres.phaseProperties(iphase, ispecies, nspecies);
        const auto cp = cres.phaseProperties(iphase, ispecies, nspecies);
        if(cp.molar_volume > 0.0)
            row(res, iphase, ispecies, nspecies) = cp.molar_volume;
        else
        {
            const auto xp = rows(x, ispecies, ispecies, nspecies, nspecies);
            row(res, iphase, ispecies, nspecies) = sum(xp % tp.standard_partial_molar_volumes);
        }

        ispecies += nspecies;
    }
    return res;
}

auto ChemicalProperties::phaseMolarEntropies() const -> VectorXd
{
    const auto& G = phaseMolarGibbsEnergies();
    const auto& H = phaseMolarEnthalpies();
    return (H - G)/T;
}

auto ChemicalProperties::phaseMolarInternalEnergies() const -> VectorXd
{
    const auto& H = phaseMolarEnthalpies();
    const auto& V = phaseMolarVolumes();
    return H - P*V;
}

auto ChemicalProperties::phaseMolarHelmholtzEnergies() const -> VectorXd
{
    const auto& G = phaseMolarGibbsEnergies();
    const auto& V = phaseMolarVolumes();
    return G - P*V;
}

auto ChemicalProperties::phaseMolarHeatCapacitiesConstP() const -> VectorXd
{
    VectorXd res(num_phases, num_species);
    Index ispecies = 0;
    for(Index iphase = 0; iphase < num_phases; ++iphase)
    {
        const auto nspecies = system.numSpeciesInPhase(iphase);
        const auto xp = rows(x, ispecies, ispecies, nspecies, nspecies);
        const auto tp = tres.phaseProperties(iphase, ispecies, nspecies);
        const auto cp = cres.phaseProperties(iphase, ispecies, nspecies);
        row(res, iphase, ispecies, nspecies) = sum(xp % tp.standard_partial_molar_heat_capacities_cp);
        row(res, iphase, ispecies, nspecies) += cp.residual_molar_heat_capacity_cp;
        ispecies += nspecies;
    }
    return res;
}

auto ChemicalProperties::phaseMolarHeatCapacitiesConstV() const -> VectorXd
{
    VectorXd res(num_phases, num_species);
    Index ispecies = 0;
    for(Index iphase = 0; iphase < num_phases; ++iphase)
    {
        const auto nspecies = system.numSpeciesInPhase(iphase);
        const auto xp = rows(x, ispecies, ispecies, nspecies, nspecies);
        const auto tp = tres.phaseProperties(iphase, ispecies, nspecies);
        const auto cp = cres.phaseProperties(iphase, ispecies, nspecies);
        row(res, iphase, ispecies, nspecies) = sum(xp % tp.standard_partial_molar_heat_capacities_cv);
        row(res, iphase, ispecies, nspecies) += cp.residual_molar_heat_capacity_cv;
        ispecies += nspecies;
    }
    return res;
}

auto ChemicalProperties::phaseSpecificGibbsEnergies() const -> VectorXd
{
    return phaseAmounts()/phaseMasses() % phaseMolarGibbsEnergies();
}

auto ChemicalProperties::phaseSpecificEnthalpies() const -> VectorXd
{
    return phaseAmounts()/phaseMasses() % phaseMolarEnthalpies();
}

auto ChemicalProperties::phaseSpecificVolumes() const -> VectorXd
{
    return phaseAmounts()/phaseMasses() % phaseMolarVolumes();
}

auto ChemicalProperties::phaseSpecificEntropies() const -> VectorXd
{
    return phaseAmounts()/phaseMasses() % phaseMolarEntropies();
}

auto ChemicalProperties::phaseSpecificInternalEnergies() const -> VectorXd
{
    return phaseAmounts()/phaseMasses() % phaseMolarInternalEnergies();
}

auto ChemicalProperties::phaseSpecificHelmholtzEnergies() const -> VectorXd
{
    return phaseAmounts()/phaseMasses() % phaseMolarHelmholtzEnergies();
}

auto ChemicalProperties::phaseSpecificHeatCapacitiesConstP() const -> VectorXd
{
    return phaseAmounts()/phaseMasses() % phaseMolarHeatCapacitiesConstP();
}

auto ChemicalProperties::phaseSpecificHeatCapacitiesConstV() const -> VectorXd
{
    return phaseAmounts()/phaseMasses() % phaseMolarHeatCapacitiesConstV();
}

auto ChemicalProperties::phaseDensities() const -> VectorXd
{
    return phaseMasses()/(phaseAmounts() % phaseMolarVolumes());
}

auto ChemicalProperties::phaseMasses() const -> VectorXd
{
    const auto mm = Reaktoro::molarMasses(system.species());
    VectorXd res(num_phases, num_species);
    Index ispecies = 0;
    for(Index iphase = 0; iphase < num_phases; ++iphase)
    {
        const auto nspecies = system.numSpeciesInPhase(iphase);
        const auto np = rows(n, ispecies, ispecies, nspecies, nspecies);
        auto mmp = rows(mm, ispecies, nspecies);
        row(res, iphase, ispecies, nspecies) = sum(mmp % np);
        ispecies += nspecies;
    }
    return res;
}

auto ChemicalProperties::phaseAmounts() const -> VectorXd
{
    VectorXd res(num_phases, num_species);
    Index ispecies = 0;
    for(Index iphase = 0; iphase < num_phases; ++iphase)
    {
        const auto nspecies = system.numSpeciesInPhase(iphase);
        const auto np = rows(n, ispecies, ispecies, nspecies, nspecies);
        row(res, iphase, ispecies, nspecies) = sum(np);
        ispecies += nspecies;
    }
    return res;
}

auto ChemicalProperties::phaseVolumes() const -> VectorXd
{
    return phaseAmounts() % phaseMolarVolumes();
}

auto ChemicalProperties::volume() const -> real
{
    return sum(phaseVolumes());
}

auto ChemicalProperties::subvolume(const Indices& iphases) const -> real
{
    return sum(rows(phaseVolumes(), iphases));
}

auto ChemicalProperties::fluidVolume() const -> real
{
    const Indices iphases = system.indicesFluidPhases();
    return sum(rows(phaseVolumes(), iphases));
}

auto ChemicalProperties::solidVolume() const -> real
{
    const Indices iphases = system.indicesSolidPhases();
    return sum(rows(phaseVolumes(), iphases));
}

} // namespace Reaktoro
