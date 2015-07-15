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

auto ChemicalProperties::temperature() const -> double
{
    return internal.T.val;
}

auto ChemicalProperties::pressure() const -> double
{
    return internal.P.val;
}

auto ChemicalProperties::composition() const -> Vector
{
    return internal.n.val;
}

auto ChemicalProperties::molarFractions() const -> ChemicalVector
{
    return internal.molar_fractions;
}

auto ChemicalProperties::lnActivityConstants() const -> ChemicalVector
{
    return internal.ln_activity_constants;
}

auto ChemicalProperties::lnActivityCoefficients() const -> ChemicalVector
{
    return internal.ln_activity_coefficients;
}

auto ChemicalProperties::lnActivities() const -> ChemicalVector
{
    return internal.ln_activities;
}

auto ChemicalProperties::chemicalPotentials() const -> ChemicalVector
{
    const auto& T = internal.T;
    const auto& R = universalGasConstant;
    const auto& G = internal.standard_partial_molar_gibbs_energies;
    const auto& lna = internal.ln_activities;
    return G + R*T*lna;
}

auto ChemicalProperties::standardPartialMolarGibbsEnergies() const -> ThermoVector
{
    return internal.standard_partial_molar_gibbs_energies;
}

auto ChemicalProperties::standardPartialMolarEnthalpies() const -> ThermoVector
{
    return internal.standard_partial_molar_enthalpies;
}

auto ChemicalProperties::standardPartialMolarVolumes() const -> ThermoVector
{
    return internal.standard_partial_molar_volumes;
}

auto ChemicalProperties::standardPartialMolarEntropies() const -> ThermoVector
{
    const auto& T = internal.T;
    const auto& G = internal.standard_partial_molar_gibbs_energies;
    const auto& H = internal.standard_partial_molar_enthalpies;
    return (H - G)/T;

}

auto ChemicalProperties::standardPartialMolarInternalEnergies() const -> ThermoVector
{
    const auto& P = internal.P;
    const auto& H = internal.standard_partial_molar_enthalpies;
    const auto& V = internal.standard_partial_molar_volumes;
    return H - P*V;
}

auto ChemicalProperties::standardPartialMolarHelmholtzEnergies() const -> ThermoVector
{
    const auto& P = internal.P;
    const auto& G = internal.standard_partial_molar_gibbs_energies;
    const auto& V = internal.standard_partial_molar_volumes;
    return G - P*V;
}

auto ChemicalProperties::standardPartialMolarHeatCapacitiesConstP() const -> ThermoVector
{
    return internal.standard_partial_molar_heat_capacities_cp;
}

auto ChemicalProperties::standardPartialMolarHeatCapacitiesConstV() const -> ThermoVector
{
    return internal.standard_partial_molar_heat_capacities_cv;
}

auto ChemicalProperties::phaseMolarGibbsEnergies() const -> ChemicalVector
{
    return internal.phase_molar_gibbs_energies;
}

auto ChemicalProperties::phaseMolarEnthalpies() const -> ChemicalVector
{
    return internal.phase_molar_enthalpies;
}

auto ChemicalProperties::phaseMolarVolumes() const -> ChemicalVector
{
    return internal.phase_molar_volumes;
}

auto ChemicalProperties::phaseMolarEntropies() const -> ChemicalVector
{
    const auto& T = internal.T;
    const auto& G = internal.phase_molar_gibbs_energies;
    const auto& H = internal.phase_molar_enthalpies;
    return (H - G)/T;
}

auto ChemicalProperties::phaseMolarInternalEnergies() const -> ChemicalVector
{
    const auto& P = internal.P;
    const auto& H = internal.phase_molar_enthalpies;
    const auto& V = internal.phase_molar_volumes;
    return H - P*V;
}

auto ChemicalProperties::phaseMolarHelmholtzEnergies() const -> ChemicalVector
{
    const auto& P = internal.P;
    const auto& G = internal.phase_molar_gibbs_energies;
    const auto& V = internal.phase_molar_volumes;
    return G - P*V;
}

auto ChemicalProperties::phaseMolarHeatCapacitiesConstP() const -> ChemicalVector
{
    return internal.phase_molar_heat_capacities_cp;
}

auto ChemicalProperties::phaseMolarHeatCapacitiesConstV() const -> ChemicalVector
{
    return internal.phase_molar_heat_capacities_cv;
}

auto ChemicalProperties::phaseSpecificGibbsEnergies() const -> ChemicalVector
{
    return internal.phase_moles/internal.phase_masses % phaseMolarGibbsEnergies();
}

auto ChemicalProperties::phaseSpecificEnthalpies() const -> ChemicalVector
{
    return internal.phase_moles/internal.phase_masses % phaseMolarEnthalpies();
}

auto ChemicalProperties::phaseSpecificVolumes() const -> ChemicalVector
{
    return internal.phase_moles/internal.phase_masses % phaseMolarVolumes();
}

auto ChemicalProperties::phaseSpecificEntropies() const -> ChemicalVector
{
    return internal.phase_moles/internal.phase_masses % phaseMolarEntropies();
}

auto ChemicalProperties::phaseSpecificInternalEnergies() const -> ChemicalVector
{
    return internal.phase_moles/internal.phase_masses % phaseMolarInternalEnergies();
}

auto ChemicalProperties::phaseSpecificHelmholtzEnergies() const -> ChemicalVector
{
    return internal.phase_moles/internal.phase_masses % phaseMolarHelmholtzEnergies();
}

auto ChemicalProperties::phaseSpecificHeatCapacitiesConstP() const -> ChemicalVector
{
    return internal.phase_moles/internal.phase_masses % phaseMolarHeatCapacitiesConstP();
}

auto ChemicalProperties::phaseSpecificHeatCapacitiesConstV() const -> ChemicalVector
{
    return internal.phase_moles/internal.phase_masses % phaseMolarHeatCapacitiesConstV();
}

auto ChemicalProperties::phaseMasses() const -> ChemicalVector
{
    return internal.phase_masses;
}

auto ChemicalProperties::phaseMoles() const -> ChemicalVector
{
    return internal.phase_moles;
}

auto ChemicalProperties::phaseVolumes() const -> ChemicalVector
{
    return internal.phase_moles % internal.phase_molar_volumes;
}


PhaseChemicalProperties::PhaseChemicalProperties()
{}

auto PhaseChemicalProperties::temperature() const -> double
{
    return internal.T.val;
}

auto PhaseChemicalProperties::pressure() const -> double
{
    return internal.P.val;
}

auto PhaseChemicalProperties::composition() const -> Vector
{
    return internal.n.val;
}

auto PhaseChemicalProperties::molarFractions() const -> ChemicalVector
{
    return internal.molar_fractions;
}

auto PhaseChemicalProperties::lnActivityConstants() const -> ChemicalVector
{
    return internal.ln_activity_constants;
}

auto PhaseChemicalProperties::lnActivityCoefficients() const -> ChemicalVector
{
    return internal.ln_activity_coefficients;
}

auto PhaseChemicalProperties::lnActivities() const -> ChemicalVector
{
    return internal.ln_activities;
}

auto PhaseChemicalProperties::chemicalPotentials() const -> ChemicalVector
{
    const auto& T = internal.T;
    const auto& R = universalGasConstant;
    const auto& G = internal.standard_partial_molar_gibbs_energies;
    const auto& lna = internal.ln_activities;
    return G + R*T*lna;
}

auto PhaseChemicalProperties::standardPartialMolarGibbsEnergies() const -> ThermoVector
{
    return internal.standard_partial_molar_gibbs_energies;
}

auto PhaseChemicalProperties::standardPartialMolarEnthalpies() const -> ThermoVector
{
    return internal.standard_partial_molar_enthalpies;
}

auto PhaseChemicalProperties::standardPartialMolarVolumes() const -> ThermoVector
{
    return internal.standard_partial_molar_volumes;
}

auto PhaseChemicalProperties::standardPartialMolarEntropies() const -> ThermoVector
{
    const auto& T = internal.T;
    const auto& G = internal.standard_partial_molar_gibbs_energies;
    const auto& H = internal.standard_partial_molar_enthalpies;
    return (H - G)/T;

}

auto PhaseChemicalProperties::standardPartialMolarInternalEnergies() const -> ThermoVector
{
    const auto& P = internal.P;
    const auto& H = internal.standard_partial_molar_enthalpies;
    const auto& V = internal.standard_partial_molar_volumes;
    return H - P*V;
}

auto PhaseChemicalProperties::standardPartialMolarHelmholtzEnergies() const -> ThermoVector
{
    const auto& P = internal.P;
    const auto& G = internal.standard_partial_molar_gibbs_energies;
    const auto& V = internal.standard_partial_molar_volumes;
    return G - P*V;
}

auto PhaseChemicalProperties::standardPartialMolarHeatCapacitiesConstP() const -> ThermoVector
{
    return internal.standard_partial_molar_heat_capacities_cp;
}

auto PhaseChemicalProperties::standardPartialMolarHeatCapacitiesConstV() const -> ThermoVector
{
    return internal.standard_partial_molar_heat_capacities_cv;
}

auto PhaseChemicalProperties::phaseMolarGibbsEnergy() const -> ChemicalScalar
{
    return internal.phase_molar_gibbs_energy;
}

auto PhaseChemicalProperties::phaseMolarEnthalpy() const -> ChemicalScalar
{
    return internal.phase_molar_enthalpy;
}

auto PhaseChemicalProperties::phaseMolarVolume() const -> ChemicalScalar
{
    return internal.phase_molar_volume;
}

auto PhaseChemicalProperties::phaseMolarEntropy() const -> ChemicalScalar
{
    const auto& T = internal.T;
    const auto& G = internal.phase_molar_gibbs_energy;
    const auto& H = internal.phase_molar_enthalpy;
    return (H - G)/T;
}

auto PhaseChemicalProperties::phaseMolarInternalEnergy() const -> ChemicalScalar
{
    const auto& P = internal.P;
    const auto& H = internal.phase_molar_enthalpy;
    const auto& V = internal.phase_molar_volume;
    return H - P*V;
}

auto PhaseChemicalProperties::phaseMolarHelmholtzEnergy() const -> ChemicalScalar
{
    const auto& P = internal.P;
    const auto& G = internal.phase_molar_gibbs_energy;
    const auto& V = internal.phase_molar_volume;
    return G - P*V;
}

auto PhaseChemicalProperties::phaseMolarHeatCapacityConstP() const -> ChemicalScalar
{
    return internal.phase_molar_heat_capacity_cp;
}

auto PhaseChemicalProperties::phaseMolarHeatCapacityConstV() const -> ChemicalScalar
{
    return internal.phase_molar_heat_capacity_cv;
}

auto PhaseChemicalProperties::phaseSpecificGibbsEnergy() const -> ChemicalScalar
{
    return internal.phase_moles/internal.phase_mass * phaseMolarGibbsEnergy();
}

auto PhaseChemicalProperties::phaseSpecificEnthalpy() const -> ChemicalScalar
{
    return internal.phase_moles/internal.phase_mass * phaseMolarEnthalpy();
}

auto PhaseChemicalProperties::phaseSpecificVolume() const -> ChemicalScalar
{
    return internal.phase_moles/internal.phase_mass * phaseMolarVolume();
}

auto PhaseChemicalProperties::phaseSpecificEntropy() const -> ChemicalScalar
{
    return internal.phase_moles/internal.phase_mass * phaseMolarEntropy();
}

auto PhaseChemicalProperties::phaseSpecificInternalEnergy() const -> ChemicalScalar
{
    return internal.phase_moles/internal.phase_mass * phaseMolarInternalEnergy();
}

auto PhaseChemicalProperties::phaseSpecificHelmholtzEnergy() const -> ChemicalScalar
{
    return internal.phase_moles/internal.phase_mass * phaseMolarHelmholtzEnergy();
}

auto PhaseChemicalProperties::phaseSpecificHeatCapacityConstP() const -> ChemicalScalar
{
    return internal.phase_moles/internal.phase_mass * phaseMolarHeatCapacityConstP();
}

auto PhaseChemicalProperties::phaseSpecificHeatCapacityConstV() const -> ChemicalScalar
{
    return internal.phase_moles/internal.phase_mass * phaseMolarHeatCapacityConstV();
}

auto PhaseChemicalProperties::phaseMass() const -> ChemicalScalar
{
    return internal.phase_mass;
}

auto PhaseChemicalProperties::phaseMoles() const -> ChemicalScalar
{
    return internal.phase_moles;
}

auto PhaseChemicalProperties::phaseVolume() const -> ChemicalScalar
{
    return internal.phase_moles * internal.phase_molar_volume;
}

} // namespace Reaktoro
