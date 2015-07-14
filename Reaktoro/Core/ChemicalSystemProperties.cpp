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

#include "ChemicalSystemProperties.hpp"

namespace Reaktoro {

ChemicalSystemProperties::ChemicalSystemProperties()
{

}

auto ChemicalSystemProperties::temperature() const -> double
{
    return T.val;
}

auto ChemicalSystemProperties::pressure() const -> double
{
    return P.val;
}

auto ChemicalSystemProperties::composition() const -> Vector
{
    return n.val;
}

auto ChemicalSystemProperties::molarFractions() const -> ChemicalVector
{
    return molar_fractions;
}

auto ChemicalSystemProperties::lnActivityConstants() const -> ChemicalVector
{
    return ln_activity_constants;
}

auto ChemicalSystemProperties::lnActivityCoefficients() const -> ChemicalVector
{
    return ln_activity_coefficients;
}

auto ChemicalSystemProperties::lnActivities() const -> ChemicalVector
{
    return ln_activities;
}

auto ChemicalSystemProperties::chemicalPotentials() const -> ChemicalVector
{
    const double R = universalGasConstant;
    return standard_partial_molar_gibbs_energies + R*T*ln_activities;
}

auto ChemicalSystemProperties::standardPartialMolarGibbsEnergies() const -> ThermoVector
{
    return standard_partial_molar_gibbs_energies;
}

auto ChemicalSystemProperties::standardPartialMolarEnthalpies() const -> ThermoVector
{
    return standard_partial_molar_enthalpies;
}

auto ChemicalSystemProperties::standardPartialMolarVolumes() const -> ThermoVector
{
    return standard_partial_molar_volumes;
}

auto ChemicalSystemProperties::standardPartialMolarEntropies() const -> ThermoVector
{
    const auto& G = standard_partial_molar_gibbs_energies;
    const auto& H = standard_partial_molar_enthalpies;
    return (H - G)/T;

}

auto ChemicalSystemProperties::standardPartialMolarInternalEnergies() const -> ThermoVector
{
    const auto& H = standard_partial_molar_enthalpies;
    const auto& V = standard_partial_molar_volumes;
    return H - P*V;
}

auto ChemicalSystemProperties::standardPartialMolarHelmholtzEnergies() const -> ThermoVector
{
    const auto& G = standard_partial_molar_gibbs_energies;
    const auto& V = standard_partial_molar_volumes;
    return G - P*V;
}

auto ChemicalSystemProperties::standardPartialMolarHeatCapacitiesConstP() const -> ThermoVector
{
    return standard_partial_molar_heat_capacities_cp;
}

auto ChemicalSystemProperties::standardPartialMolarHeatCapacitiesConstV() const -> ThermoVector
{
    return standard_partial_molar_heat_capacities_cv;
}

auto ChemicalSystemProperties::phaseMolarGibbsEnergies() const -> ChemicalVector
{
    return phase_molar_gibbs_energies;
}

auto ChemicalSystemProperties::phaseMolarEnthalpies() const -> ChemicalVector
{
    return phase_molar_enthalpies;
}

auto ChemicalSystemProperties::phaseMolarVolumes() const -> ChemicalVector
{
    return phase_molar_volumes;
}

auto ChemicalSystemProperties::phaseMolarEntropies() const -> ChemicalVector
{
    const auto& G = phase_molar_gibbs_energies;
    const auto& H = phase_molar_enthalpies;
    return (H - G)/T;
}

auto ChemicalSystemProperties::phaseMolarInternalEnergies() const -> ChemicalVector
{
    const auto& H = phase_molar_enthalpies;
    const auto& V = phase_molar_volumes;
    return H - P*V;
}

auto ChemicalSystemProperties::phaseMolarHelmholtzEnergies() const -> ChemicalVector
{
    const auto& G = phase_molar_gibbs_energies;
    const auto& V = phase_molar_volumes;
    return G - P*V;
}

auto ChemicalSystemProperties::phaseMolarHeatCapacitiesConstP() const -> ChemicalVector
{
    return phase_molar_heat_capacities_cp;
}

auto ChemicalSystemProperties::phaseMolarHeatCapacitiesConstV() const -> ChemicalVector
{
    return phase_molar_heat_capacities_cv;
}

auto ChemicalSystemProperties::phaseSpecificGibbsEnergies() const -> ChemicalVector
{
    return phase_moles/phase_masses % phaseMolarGibbsEnergies();
}

auto ChemicalSystemProperties::phaseSpecificEnthalpies() const -> ChemicalVector
{
    return phase_moles/phase_masses % phaseMolarEnthalpies();
}

auto ChemicalSystemProperties::phaseSpecificVolumes() const -> ChemicalVector
{
    return phase_moles/phase_masses % phaseMolarVolumes();
}

auto ChemicalSystemProperties::phaseSpecificEntropies() const -> ChemicalVector
{
    return phase_moles/phase_masses % phaseMolarEntropies();
}

auto ChemicalSystemProperties::phaseSpecificInternalEnergies() const -> ChemicalVector
{
    return phase_moles/phase_masses % phaseMolarInternalEnergies();
}

auto ChemicalSystemProperties::phaseSpecificHelmholtzEnergies() const -> ChemicalVector
{
    return phase_moles/phase_masses % phaseMolarHelmholtzEnergies();
}

auto ChemicalSystemProperties::phaseSpecificHeatCapacitiesConstP() const -> ChemicalVector
{
    return phase_moles/phase_masses % phaseMolarHeatCapacitiesConstP();
}

auto ChemicalSystemProperties::phaseSpecificHeatCapacitiesConstV() const -> ChemicalVector
{
    return phase_moles/phase_masses % phaseMolarHeatCapacitiesConstV();
}

auto ChemicalSystemProperties::phaseMasses() const -> ChemicalVector
{
    return phase_masses;
}

auto ChemicalSystemProperties::phaseMoles() const -> ChemicalVector
{
    return phase_moles;
}

auto ChemicalSystemProperties::phaseVolumes() const -> ChemicalVector
{
    return phase_moles % phase_molar_volumes;
}

} // namespace Reaktoro
