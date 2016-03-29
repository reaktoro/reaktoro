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
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {

ChemicalProperties::ChemicalProperties()
{}

ChemicalProperties::ChemicalProperties(const ChemicalSystem& system)
: system(system)
{}

auto ChemicalProperties::temperature() const -> double
{
    return T.val;
}

auto ChemicalProperties::pressure() const -> double
{
    return P.val;
}

auto ChemicalProperties::composition() const -> const Vector&
{
    return n;
}

auto ChemicalProperties::molarFractions() const -> ChemicalVector
{
    return x;
}

auto ChemicalProperties::lnActivityCoefficients() const -> ChemicalVector
{
    ChemicalVector res(x.size());
    unsigned offset = 0;
    for(unsigned i = 0; i < cresults.size(); ++i)
    {
        const unsigned size = cresults[i].num_species;
        res.rows(offset, offset, size, size) = cresults[i].ln_activity_coefficients;
        offset += size;
    }
    return res;
}

auto ChemicalProperties::lnActivityConstants() const -> ThermoVector
{
    ThermoVector res(x.size());
    unsigned offset = 0;
    for(unsigned i = 0; i < cresults.size(); ++i)
    {
        const unsigned size = cresults[i].num_species;
        res.rows(offset, size) = cresults[i].ln_activity_constants;
        offset += size;
    }
    return res;
}

auto ChemicalProperties::lnActivities() const -> ChemicalVector
{
    ChemicalVector res(x.size());
    unsigned offset = 0;
    for(unsigned i = 0; i < cresults.size(); ++i)
    {
        const unsigned size = cresults[i].num_species;
        res.rows(offset, offset, size, size) = cresults[i].ln_activities;
        offset += size;
    }
    return res;
}

auto ChemicalProperties::chemicalPotentials() const -> ChemicalVector
{
    const auto& R = universalGasConstant;
    const auto& G = standardPartialMolarGibbsEnergies();
    const auto& lna = lnActivities();
    return G + R*T*lna;
}

auto ChemicalProperties::standardPartialMolarGibbsEnergies() const -> ThermoVector
{
    ThermoVector res(x.size());
    unsigned offset = 0;
    for(unsigned i = 0; i < tresults.size(); ++i)
    {
        const unsigned size = tresults[i].num_species;
        res.rows(offset, size) = tresults[i].standard_partial_molar_gibbs_energies;
        offset += size;
    }
    return res;
}

auto ChemicalProperties::standardPartialMolarEnthalpies() const -> ThermoVector
{
    ThermoVector res(x.size());
    unsigned offset = 0;
    for(unsigned i = 0; i < tresults.size(); ++i)
    {
        const unsigned size = tresults[i].num_species;
        res.rows(offset, size) = tresults[i].standard_partial_molar_enthalpies;
        offset += size;
    }
    return res;
}

auto ChemicalProperties::standardPartialMolarVolumes() const -> ThermoVector
{
    ThermoVector res(x.size());
    unsigned offset = 0;
    for(unsigned i = 0; i < tresults.size(); ++i)
    {
        const unsigned size = tresults[i].num_species;
        res.rows(offset, size) = tresults[i].standard_partial_molar_volumes;
        offset += size;
    }
    return res;
}

auto ChemicalProperties::standardPartialMolarEntropies() const -> ThermoVector
{
    const auto& G = standardPartialMolarGibbsEnergies();
    const auto& H = standardPartialMolarEnthalpies();
    return (H - G)/T;

}

auto ChemicalProperties::standardPartialMolarInternalEnergies() const -> ThermoVector
{
    const auto& H = standardPartialMolarEnthalpies();
    const auto& V = standardPartialMolarVolumes();
    return H - P*V;
}

auto ChemicalProperties::standardPartialMolarHelmholtzEnergies() const -> ThermoVector
{
    const auto& G = standardPartialMolarGibbsEnergies();
    const auto& V = standardPartialMolarVolumes();
    return G - P*V;
}

auto ChemicalProperties::standardPartialMolarHeatCapacitiesConstP() const -> ThermoVector
{
    ThermoVector res(x.size());
    unsigned offset = 0;
    for(unsigned i = 0; i < tresults.size(); ++i)
    {
        const unsigned size = tresults[i].num_species;
        res.rows(offset, size) = tresults[i].standard_partial_molar_heat_capacities_cp;
        offset += size;
    }
    return res;
}

auto ChemicalProperties::standardPartialMolarHeatCapacitiesConstV() const -> ThermoVector
{
    ThermoVector res(x.size());
    unsigned offset = 0;
    for(unsigned i = 0; i < tresults.size(); ++i)
    {
        const unsigned size = tresults[i].num_species;
        res.rows(offset, size) = tresults[i].standard_partial_molar_heat_capacities_cv;
        offset += size;
    }
    return res;
}

auto ChemicalProperties::phaseMolarGibbsEnergies() const -> ChemicalVector
{
    ChemicalVector res(cresults.size(), x.size());
    unsigned offset = 0;
    for(unsigned i = 0; i < cresults.size(); ++i)
    {
        const unsigned size = cresults[i].num_species;
        const auto xp = x.rows(offset, offset, size, size);
        res.row(i, offset, size) = sum(xp % tresults[i].standard_partial_molar_gibbs_energies);
        offset += size;
    }
    return res;
}

auto ChemicalProperties::phaseMolarEnthalpies() const -> ChemicalVector
{
    ChemicalVector res(cresults.size(), x.size());
    unsigned offset = 0;
    for(unsigned i = 0; i < cresults.size(); ++i)
    {
        const unsigned size = cresults[i].num_species;
        const auto xp = x.rows(offset, offset, size, size);
        res.row(i, offset, size) = sum(xp % tresults[i].standard_partial_molar_enthalpies);
        offset += size;
    }
    return res;
}

auto ChemicalProperties::phaseMolarVolumes() const -> ChemicalVector
{
    ChemicalVector res(cresults.size(), x.size());
    unsigned offset = 0;
    for(unsigned i = 0; i < cresults.size(); ++i)
    {
        const unsigned size = cresults[i].num_species;
        const auto xp = x.rows(offset, offset, size, size);
        if(cresults[i].molar_volume.val > 0.0)
            res.row(i, offset, size) = cresults[i].molar_volume;
        else
            res.row(i, offset, size) = sum(xp % tresults[i].standard_partial_molar_volumes);

        offset += size;
    }
    return res;
}

auto ChemicalProperties::phaseMolarEntropies() const -> ChemicalVector
{
    const auto& G = phaseMolarGibbsEnergies();
    const auto& H = phaseMolarEnthalpies();
    return (H - G)/T;
}

auto ChemicalProperties::phaseMolarInternalEnergies() const -> ChemicalVector
{
    const auto& H = phaseMolarEnthalpies();
    const auto& V = phaseMolarVolumes();
    return H - P*V;
}

auto ChemicalProperties::phaseMolarHelmholtzEnergies() const -> ChemicalVector
{
    const auto& G = phaseMolarGibbsEnergies();
    const auto& V = phaseMolarVolumes();
    return G - P*V;
}

auto ChemicalProperties::phaseMolarHeatCapacitiesConstP() const -> ChemicalVector
{
    ChemicalVector res(cresults.size(), x.size());
    unsigned offset = 0;
    for(unsigned i = 0; i < cresults.size(); ++i)
    {
        const unsigned size = cresults[i].num_species;
        const auto xp = x.rows(offset, offset, size, size);
        res.row(i, offset, size) = sum(xp % tresults[i].standard_partial_molar_heat_capacities_cp);
        offset += size;
    }
    return res;
}

auto ChemicalProperties::phaseMolarHeatCapacitiesConstV() const -> ChemicalVector
{
    ChemicalVector res(cresults.size(), x.size());
    unsigned offset = 0;
    for(unsigned i = 0; i < cresults.size(); ++i)
    {
        const unsigned size = cresults[i].num_species;
        const auto xp = x.rows(offset, offset, size, size);
        res.row(i, offset, size) = sum(xp % tresults[i].standard_partial_molar_heat_capacities_cv);
        offset += size;
    }
    return res;
}

auto ChemicalProperties::phaseSpecificGibbsEnergies() const -> ChemicalVector
{
    return phaseAmounts()/phaseMasses() % phaseMolarGibbsEnergies();
}

auto ChemicalProperties::phaseSpecificEnthalpies() const -> ChemicalVector
{
    return phaseAmounts()/phaseMasses() % phaseMolarEnthalpies();
}

auto ChemicalProperties::phaseSpecificVolumes() const -> ChemicalVector
{
    return phaseAmounts()/phaseMasses() % phaseMolarVolumes();
}

auto ChemicalProperties::phaseSpecificEntropies() const -> ChemicalVector
{
    return phaseAmounts()/phaseMasses() % phaseMolarEntropies();
}

auto ChemicalProperties::phaseSpecificInternalEnergies() const -> ChemicalVector
{
    return phaseAmounts()/phaseMasses() % phaseMolarInternalEnergies();
}

auto ChemicalProperties::phaseSpecificHelmholtzEnergies() const -> ChemicalVector
{
    return phaseAmounts()/phaseMasses() % phaseMolarHelmholtzEnergies();
}

auto ChemicalProperties::phaseSpecificHeatCapacitiesConstP() const -> ChemicalVector
{
    return phaseAmounts()/phaseMasses() % phaseMolarHeatCapacitiesConstP();
}

auto ChemicalProperties::phaseSpecificHeatCapacitiesConstV() const -> ChemicalVector
{
    return phaseAmounts()/phaseMasses() % phaseMolarHeatCapacitiesConstV();
}

auto ChemicalProperties::phaseDensities() const -> ChemicalVector
{
    return phaseMasses()/(phaseAmounts() % phaseMolarVolumes());
}

auto ChemicalProperties::phaseMasses() const -> ChemicalVector
{
    auto nc = Reaktoro::composition(n);
    auto mm = Reaktoro::molarMasses(system.species());
    ChemicalVector res(cresults.size(), x.size());
    unsigned offset = 0;
    for(unsigned i = 0; i < cresults.size(); ++i)
    {
        const unsigned size = cresults[i].num_species;
        auto np = nc.rows(offset, offset, size, size);
        auto mmp = rows(mm, offset, size);
        res.row(i, offset, size) = sum(mmp % np);
        offset += size;
    }
    return res;
}

auto ChemicalProperties::phaseAmounts() const -> ChemicalVector
{
    auto nc = Reaktoro::composition(n);
    ChemicalVector res(cresults.size(), x.size());
    unsigned offset = 0;
    for(unsigned i = 0; i < cresults.size(); ++i)
    {
        const unsigned size = cresults[i].num_species;
        auto np = nc.rows(offset, offset, size, size);
        res.row(i, offset, size) = sum(np);
        offset += size;
    }
    return res;
}

auto ChemicalProperties::phaseVolumes() const -> ChemicalVector
{
    return phaseAmounts() % phaseMolarVolumes();
}

auto ChemicalProperties::volume() const -> ChemicalScalar
{
    return sum(phaseAmounts() % phaseMolarVolumes());
}

PhaseChemicalProperties::PhaseChemicalProperties()
{}

PhaseChemicalProperties::PhaseChemicalProperties(unsigned nspecies)
: standard_partial_molar_gibbs_energies(nspecies),
  standard_partial_molar_enthalpies(nspecies),
  standard_partial_molar_volumes(nspecies),
  standard_partial_molar_heat_capacities_cp(nspecies),
  standard_partial_molar_heat_capacities_cv(nspecies),
  molar_fractions(nspecies),
  ln_activity_coefficients(nspecies),
  ln_activity_constants(nspecies),
  ln_activities(nspecies),
  phase_molar_gibbs_energy(nspecies),
  phase_molar_enthalpy(nspecies),
  phase_molar_volume(nspecies),
  phase_molar_heat_capacity_cp(nspecies),
  phase_molar_heat_capacity_cv(nspecies),
  phase_amount(nspecies),
  phase_mass(nspecies)
{}

auto PhaseChemicalProperties::temperature() const -> double
{
    return T.val;
}

auto PhaseChemicalProperties::pressure() const -> double
{
    return P.val;
}

auto PhaseChemicalProperties::composition() const -> Vector
{
    return n;
}

auto PhaseChemicalProperties::molarFractions() const -> ChemicalVector
{
    return molar_fractions;
}

auto PhaseChemicalProperties::lnActivityCoefficients() const -> ChemicalVector
{
    return ln_activity_coefficients;
}

auto PhaseChemicalProperties::lnActivityConstants() const -> ThermoVector
{
    return ln_activity_constants;
}

auto PhaseChemicalProperties::lnActivities() const -> ChemicalVector
{
    return ln_activities;
}

auto PhaseChemicalProperties::chemicalPotentials() const -> ChemicalVector
{
    const auto& R = universalGasConstant;
    const auto& G = standard_partial_molar_gibbs_energies;
    const auto& lna = ln_activities;
    return G + R*T*lna;
}

auto PhaseChemicalProperties::standardPartialMolarGibbsEnergies() const -> ThermoVector
{
    return standard_partial_molar_gibbs_energies;
}

auto PhaseChemicalProperties::standardPartialMolarEnthalpies() const -> ThermoVector
{
    return standard_partial_molar_enthalpies;
}

auto PhaseChemicalProperties::standardPartialMolarVolumes() const -> ThermoVector
{
    return standard_partial_molar_volumes;
}

auto PhaseChemicalProperties::standardPartialMolarEntropies() const -> ThermoVector
{
    const auto& G = standard_partial_molar_gibbs_energies;
    const auto& H = standard_partial_molar_enthalpies;
    return (H - G)/T;

}

auto PhaseChemicalProperties::standardPartialMolarInternalEnergies() const -> ThermoVector
{
    const auto& H = standard_partial_molar_enthalpies;
    const auto& V = standard_partial_molar_volumes;
    return H - P*V;
}

auto PhaseChemicalProperties::standardPartialMolarHelmholtzEnergies() const -> ThermoVector
{
    const auto& G = standard_partial_molar_gibbs_energies;
    const auto& V = standard_partial_molar_volumes;
    return G - P*V;
}

auto PhaseChemicalProperties::standardPartialMolarHeatCapacitiesConstP() const -> ThermoVector
{
    return standard_partial_molar_heat_capacities_cp;
}

auto PhaseChemicalProperties::standardPartialMolarHeatCapacitiesConstV() const -> ThermoVector
{
    return standard_partial_molar_heat_capacities_cv;
}

auto PhaseChemicalProperties::molarGibbsEnergy() const -> ChemicalScalar
{
    return phase_molar_gibbs_energy;
}

auto PhaseChemicalProperties::molarEnthalpy() const -> ChemicalScalar
{
    return phase_molar_enthalpy;
}

auto PhaseChemicalProperties::molarVolume() const -> ChemicalScalar
{
    return phase_molar_volume;
}

auto PhaseChemicalProperties::molarEntropy() const -> ChemicalScalar
{
    const auto& G = phase_molar_gibbs_energy;
    const auto& H = phase_molar_enthalpy;
    return (H - G)/T;
}

auto PhaseChemicalProperties::molarInternalEnergy() const -> ChemicalScalar
{
    const auto& H = phase_molar_enthalpy;
    const auto& V = phase_molar_volume;
    return H - P*V;
}

auto PhaseChemicalProperties::molarHelmholtzEnergy() const -> ChemicalScalar
{
    const auto& G = phase_molar_gibbs_energy;
    const auto& V = phase_molar_volume;
    return G - P*V;
}

auto PhaseChemicalProperties::molarHeatCapacityConstP() const -> ChemicalScalar
{
    return phase_molar_heat_capacity_cp;
}

auto PhaseChemicalProperties::molarHeatCapacityConstV() const -> ChemicalScalar
{
    return phase_molar_heat_capacity_cv;
}

auto PhaseChemicalProperties::specificGibbsEnergy() const -> ChemicalScalar
{
    return phase_amount/phase_mass * molarGibbsEnergy();
}

auto PhaseChemicalProperties::specificEnthalpy() const -> ChemicalScalar
{
    return phase_amount/phase_mass * molarEnthalpy();
}

auto PhaseChemicalProperties::specificVolume() const -> ChemicalScalar
{
    return phase_amount/phase_mass * molarVolume();
}

auto PhaseChemicalProperties::specificEntropy() const -> ChemicalScalar
{
    return phase_amount/phase_mass * molarEntropy();
}

auto PhaseChemicalProperties::specificInternalEnergy() const -> ChemicalScalar
{
    return phase_amount/phase_mass * molarInternalEnergy();
}

auto PhaseChemicalProperties::specificHelmholtzEnergy() const -> ChemicalScalar
{
    return phase_amount/phase_mass * molarHelmholtzEnergy();
}

auto PhaseChemicalProperties::specificHeatCapacityConstP() const -> ChemicalScalar
{
    return phase_amount/phase_mass * molarHeatCapacityConstP();
}

auto PhaseChemicalProperties::specificHeatCapacityConstV() const -> ChemicalScalar
{
    return phase_amount/phase_mass * molarHeatCapacityConstV();
}

auto PhaseChemicalProperties::density() const -> ChemicalScalar
{
    return phase_mass/(phase_amount * phase_molar_volume);
}

auto PhaseChemicalProperties::mass() const -> ChemicalScalar
{
    return phase_mass;
}

auto PhaseChemicalProperties::amount() const -> ChemicalScalar
{
    return phase_amount;
}

auto PhaseChemicalProperties::volume() const -> ChemicalScalar
{
    return phase_amount * phase_molar_volume;
}

} // namespace Reaktoro
