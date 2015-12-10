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

#include "ChemicalProperties.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>

namespace Reaktoro {

ChemicalProperties::ChemicalProperties()
{}

ChemicalProperties::ChemicalProperties(unsigned nspecies, unsigned nphases)
: standard_partial_molar_gibbs_energies(nspecies),
  standard_partial_molar_enthalpies(nspecies),
  standard_partial_molar_volumes(nspecies),
  standard_partial_molar_heat_capacities_cp(nspecies),
  standard_partial_molar_heat_capacities_cv(nspecies),
  molar_fractions(nspecies),
  ln_activity_coefficients(nspecies),
  ln_activities(nspecies),
  phase_molar_gibbs_energies(nphases, nspecies),
  phase_molar_enthalpies(nphases, nspecies),
  phase_molar_volumes(nphases, nspecies),
  phase_molar_heat_capacities_cp(nphases, nspecies),
  phase_molar_heat_capacities_cv(nphases, nspecies),
  phase_moles(nphases, nspecies),
  phase_masses(nphases, nspecies)
{}

auto ChemicalProperties::temperature() const -> double
{
    return T.val;
}

auto ChemicalProperties::pressure() const -> double
{
    return P.val;
}

auto ChemicalProperties::composition() const -> Vector
{
    return n;
}

auto ChemicalProperties::molarFractions() const -> ChemicalVector
{
    return molar_fractions;
}

auto ChemicalProperties::lnActivityCoefficients() const -> ChemicalVector
{
    return ln_activity_coefficients;
}

auto ChemicalProperties::lnActivities() const -> ChemicalVector
{
    return ln_activities;
}

auto ChemicalProperties::chemicalPotentials() const -> ChemicalVector
{
    const auto& R = universalGasConstant;
    const auto& G = standard_partial_molar_gibbs_energies;
    const auto& lna = ln_activities;
    return G + R*T*lna;
}

auto ChemicalProperties::standardPartialMolarGibbsEnergies() const -> ThermoVector
{
    return standard_partial_molar_gibbs_energies;
}

auto ChemicalProperties::standardPartialMolarEnthalpies() const -> ThermoVector
{
    return standard_partial_molar_enthalpies;
}

auto ChemicalProperties::standardPartialMolarVolumes() const -> ThermoVector
{
    return standard_partial_molar_volumes;
}

auto ChemicalProperties::standardPartialMolarEntropies() const -> ThermoVector
{
    const auto& G = standard_partial_molar_gibbs_energies;
    const auto& H = standard_partial_molar_enthalpies;
    return (H - G)/T;

}

auto ChemicalProperties::standardPartialMolarInternalEnergies() const -> ThermoVector
{
    const auto& H = standard_partial_molar_enthalpies;
    const auto& V = standard_partial_molar_volumes;
    return H - P*V;
}

auto ChemicalProperties::standardPartialMolarHelmholtzEnergies() const -> ThermoVector
{
    const auto& G = standard_partial_molar_gibbs_energies;
    const auto& V = standard_partial_molar_volumes;
    return G - P*V;
}

auto ChemicalProperties::standardPartialMolarHeatCapacitiesConstP() const -> ThermoVector
{
    return standard_partial_molar_heat_capacities_cp;
}

auto ChemicalProperties::standardPartialMolarHeatCapacitiesConstV() const -> ThermoVector
{
    return standard_partial_molar_heat_capacities_cv;
}

auto ChemicalProperties::phaseMolarGibbsEnergies() const -> ChemicalVector
{
    return phase_molar_gibbs_energies;
}

auto ChemicalProperties::phaseMolarEnthalpies() const -> ChemicalVector
{
    return phase_molar_enthalpies;
}

auto ChemicalProperties::phaseMolarVolumes() const -> ChemicalVector
{
    return phase_molar_volumes;
}

auto ChemicalProperties::phaseMolarEntropies() const -> ChemicalVector
{
    const auto& G = phase_molar_gibbs_energies;
    const auto& H = phase_molar_enthalpies;
    return (H - G)/T;
}

auto ChemicalProperties::phaseMolarInternalEnergies() const -> ChemicalVector
{
    const auto& H = phase_molar_enthalpies;
    const auto& V = phase_molar_volumes;
    return H - P*V;
}

auto ChemicalProperties::phaseMolarHelmholtzEnergies() const -> ChemicalVector
{
    const auto& G = phase_molar_gibbs_energies;
    const auto& V = phase_molar_volumes;
    return G - P*V;
}

auto ChemicalProperties::phaseMolarHeatCapacitiesConstP() const -> ChemicalVector
{
    return phase_molar_heat_capacities_cp;
}

auto ChemicalProperties::phaseMolarHeatCapacitiesConstV() const -> ChemicalVector
{
    return phase_molar_heat_capacities_cv;
}

auto ChemicalProperties::phaseSpecificGibbsEnergies() const -> ChemicalVector
{
    return phase_moles/phase_masses % phaseMolarGibbsEnergies();
}

auto ChemicalProperties::phaseSpecificEnthalpies() const -> ChemicalVector
{
    return phase_moles/phase_masses % phaseMolarEnthalpies();
}

auto ChemicalProperties::phaseSpecificVolumes() const -> ChemicalVector
{
    return phase_moles/phase_masses % phaseMolarVolumes();
}

auto ChemicalProperties::phaseSpecificEntropies() const -> ChemicalVector
{
    return phase_moles/phase_masses % phaseMolarEntropies();
}

auto ChemicalProperties::phaseSpecificInternalEnergies() const -> ChemicalVector
{
    return phase_moles/phase_masses % phaseMolarInternalEnergies();
}

auto ChemicalProperties::phaseSpecificHelmholtzEnergies() const -> ChemicalVector
{
    return phase_moles/phase_masses % phaseMolarHelmholtzEnergies();
}

auto ChemicalProperties::phaseSpecificHeatCapacitiesConstP() const -> ChemicalVector
{
    return phase_moles/phase_masses % phaseMolarHeatCapacitiesConstP();
}

auto ChemicalProperties::phaseSpecificHeatCapacitiesConstV() const -> ChemicalVector
{
    return phase_moles/phase_masses % phaseMolarHeatCapacitiesConstV();
}

auto ChemicalProperties::phaseMasses() const -> ChemicalVector
{
    return phase_masses;
}

auto ChemicalProperties::phaseMoles() const -> ChemicalVector
{
    return phase_moles;
}

auto ChemicalProperties::phaseVolumes() const -> ChemicalVector
{
    return phase_moles % phase_molar_volumes;
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
  ln_activities(nspecies),
  phase_molar_gibbs_energy(nspecies),
  phase_molar_enthalpy(nspecies),
  phase_molar_volume(nspecies),
  phase_molar_heat_capacity_cp(nspecies),
  phase_molar_heat_capacity_cv(nspecies),
  phase_moles(nspecies),
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

auto PhaseChemicalProperties::phaseMolarGibbsEnergy() const -> ChemicalScalar
{
    return phase_molar_gibbs_energy;
}

auto PhaseChemicalProperties::phaseMolarEnthalpy() const -> ChemicalScalar
{
    return phase_molar_enthalpy;
}

auto PhaseChemicalProperties::phaseMolarVolume() const -> ChemicalScalar
{
    return phase_molar_volume;
}

auto PhaseChemicalProperties::phaseMolarEntropy() const -> ChemicalScalar
{
    const auto& G = phase_molar_gibbs_energy;
    const auto& H = phase_molar_enthalpy;
    return (H - G)/T;
}

auto PhaseChemicalProperties::phaseMolarInternalEnergy() const -> ChemicalScalar
{
    const auto& H = phase_molar_enthalpy;
    const auto& V = phase_molar_volume;
    return H - P*V;
}

auto PhaseChemicalProperties::phaseMolarHelmholtzEnergy() const -> ChemicalScalar
{
    const auto& G = phase_molar_gibbs_energy;
    const auto& V = phase_molar_volume;
    return G - P*V;
}

auto PhaseChemicalProperties::phaseMolarHeatCapacityConstP() const -> ChemicalScalar
{
    return phase_molar_heat_capacity_cp;
}

auto PhaseChemicalProperties::phaseMolarHeatCapacityConstV() const -> ChemicalScalar
{
    return phase_molar_heat_capacity_cv;
}

auto PhaseChemicalProperties::phaseSpecificGibbsEnergy() const -> ChemicalScalar
{
    return phase_moles/phase_mass * phaseMolarGibbsEnergy();
}

auto PhaseChemicalProperties::phaseSpecificEnthalpy() const -> ChemicalScalar
{
    return phase_moles/phase_mass * phaseMolarEnthalpy();
}

auto PhaseChemicalProperties::phaseSpecificVolume() const -> ChemicalScalar
{
    return phase_moles/phase_mass * phaseMolarVolume();
}

auto PhaseChemicalProperties::phaseSpecificEntropy() const -> ChemicalScalar
{
    return phase_moles/phase_mass * phaseMolarEntropy();
}

auto PhaseChemicalProperties::phaseSpecificInternalEnergy() const -> ChemicalScalar
{
    return phase_moles/phase_mass * phaseMolarInternalEnergy();
}

auto PhaseChemicalProperties::phaseSpecificHelmholtzEnergy() const -> ChemicalScalar
{
    return phase_moles/phase_mass * phaseMolarHelmholtzEnergy();
}

auto PhaseChemicalProperties::phaseSpecificHeatCapacityConstP() const -> ChemicalScalar
{
    return phase_moles/phase_mass * phaseMolarHeatCapacityConstP();
}

auto PhaseChemicalProperties::phaseSpecificHeatCapacityConstV() const -> ChemicalScalar
{
    return phase_moles/phase_mass * phaseMolarHeatCapacityConstV();
}

auto PhaseChemicalProperties::phaseMass() const -> ChemicalScalar
{
    return phase_mass;
}

auto PhaseChemicalProperties::phaseMoles() const -> ChemicalScalar
{
    return phase_moles;
}

auto PhaseChemicalProperties::phaseVolume() const -> ChemicalScalar
{
    return phase_moles * phase_molar_volume;
}

} // namespace Reaktoro
