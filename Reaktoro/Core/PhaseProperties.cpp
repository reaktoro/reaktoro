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

#include "PhaseProperties.hpp"

namespace Reaktoro {

PhaseProperties::PhaseProperties()
{

}

auto PhaseProperties::temperature() const -> double
{
    return T.val;
}

auto PhaseProperties::pressure() const -> double
{
    return P.val;
}

auto PhaseProperties::composition() const -> Vector
{
    return n.val;
}

auto PhaseProperties::molarFractions() const -> ChemicalVector
{
    return molar_fractions;
}

auto PhaseProperties::lnActivityConstants() const -> ChemicalVector
{
    return ln_activity_constants;
}

auto PhaseProperties::lnActivityCoefficients() const -> ChemicalVector
{
    return ln_activity_coefficients;
}

auto PhaseProperties::lnActivities() const -> ChemicalVector
{
    return ln_activities;
}

auto PhaseProperties::standardPartialMolarGibbsEnergies() const -> ThermoVector
{
    return standard_partial_molar_gibbs_energies;
}

auto PhaseProperties::standardPartialMolarEnthalpies() const -> ThermoVector
{
    return standard_partial_molar_enthalpies;
}

auto PhaseProperties::standardPartialMolarVolumes() const -> ThermoVector
{
    return standard_partial_molar_volumes;
}

auto PhaseProperties::standardPartialMolarEntropies() const -> ThermoVector
{
    const auto& G = standard_partial_molar_gibbs_energies;
    const auto& H = standard_partial_molar_enthalpies;
    return (H - G)/T;

}

auto PhaseProperties::standardPartialMolarInternalEnergies() const -> ThermoVector
{
    const auto& H = standard_partial_molar_enthalpies;
    const auto& V = standard_partial_molar_volumes;
    return H - P*V;
}

auto PhaseProperties::standardPartialMolarHelmholtzEnergies() const -> ThermoVector
{
    const auto& G = standard_partial_molar_gibbs_energies;
    const auto& V = standard_partial_molar_volumes;
    return G - P*V;
}

auto PhaseProperties::standardPartialMolarHeatCapacitiesConstP() const -> ThermoVector
{
    return standard_partial_molar_heat_capacities_cp;
}

auto PhaseProperties::standardPartialMolarHeatCapacitiesConstV() const -> ThermoVector
{
    return standard_partial_molar_heat_capacities_cv;
}

auto PhaseProperties::molarGibbsEnergy() const -> ChemicalScalar
{
    return molar_gibbs_energy;
}

auto PhaseProperties::molarEnthalpy() const -> ChemicalScalar
{
    return molar_enthalpy;
}

auto PhaseProperties::molarVolume() const -> ChemicalScalar
{
    return molar_volume;
}

auto PhaseProperties::molarEntropy() const -> ChemicalScalar
{
    const auto& G = molar_gibbs_energy;
    const auto& H = molar_enthalpy;
    return (H - G)/T;
}

auto PhaseProperties::molarInternalEnergy() const -> ChemicalScalar
{
    const auto& H = molar_enthalpy;
    const auto& V = molar_volume;
    return H - P*V;
}

auto PhaseProperties::molarHelmholtzEnergy() const -> ChemicalScalar
{
    const auto& G = molar_gibbs_energy;
    const auto& V = molar_volume;
    return G - P*V;
}

auto PhaseProperties::molarHeatCapacityConstP() const -> ChemicalScalar
{
    return molar_heat_capacity_cp;
}

auto PhaseProperties::molarHeatCapacityConstV() const -> ChemicalScalar
{
    return molar_heat_capacity_cv;
}

auto PhaseProperties::specificGibbsEnergy() const -> ChemicalScalar
{
    return moles()/mass() * molarGibbsEnergy();
}

auto PhaseProperties::specificEnthalpy() const -> ChemicalScalar
{
    return moles()/mass() * molarEnthalpy();
}

auto PhaseProperties::specificVolume() const -> ChemicalScalar
{
    return moles()/mass() * molarVolume();
}

auto PhaseProperties::specificEntropy() const -> ChemicalScalar
{
    return moles()/mass() * molarEntropy();
}

auto PhaseProperties::specificInternalEnergy() const -> ChemicalScalar
{
    return moles()/mass() * molarInternalEnergy();
}

auto PhaseProperties::specificHelmholtzEnergy() const -> ChemicalScalar
{
    return moles()/mass() * molarHelmholtzEnergy();
}

auto PhaseProperties::specificHeatCapacityConstP() const -> ChemicalScalar
{
    return moles()/mass() * molarHeatCapacityConstP();
}

auto PhaseProperties::specificHeatCapacityConstV() const -> ChemicalScalar
{
    return moles()/mass() * molarHeatCapacityConstV();
}

auto PhaseProperties::moles() const -> ChemicalScalar
{
    return sum(n);
}

auto PhaseProperties::mass() const -> ChemicalScalar
{
    return total_mass;
}

} // namespace Reaktoro
