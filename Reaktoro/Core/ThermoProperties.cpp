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

#include "ThermoProperties.hpp"

namespace Reaktoro {

ThermoProperties::ThermoProperties()
{}

ThermoProperties::ThermoProperties(unsigned nspecies)
: standard_partial_molar_gibbs_energies(nspecies),
  standard_partial_molar_enthalpies(nspecies),
  standard_partial_molar_volumes(nspecies),
  standard_partial_molar_heat_capacities_cp(nspecies),
  standard_partial_molar_heat_capacities_cv(nspecies)
{}

auto ThermoProperties::temperature() const -> double
{
    return T;
}

auto ThermoProperties::pressure() const -> double
{
    return P;
}

auto ThermoProperties::standardPartialMolarGibbsEnergies() const -> ThermoVector
{
    return standard_partial_molar_gibbs_energies;
}

auto ThermoProperties::standardPartialMolarEnthalpies() const -> ThermoVector
{
    return standard_partial_molar_enthalpies;
}

auto ThermoProperties::standardPartialMolarVolumes() const -> ThermoVector
{
    return standard_partial_molar_volumes;
}

auto ThermoProperties::standardPartialMolarEntropies() const -> ThermoVector
{
    const auto& G = standard_partial_molar_gibbs_energies;
    const auto& H = standard_partial_molar_enthalpies;
    return (H - G)/T;

}

auto ThermoProperties::standardPartialMolarInternalEnergies() const -> ThermoVector
{
    const auto& H = standard_partial_molar_enthalpies;
    const auto& V = standard_partial_molar_volumes;
    return H - P*V;
}

auto ThermoProperties::standardPartialMolarHelmholtzEnergies() const -> ThermoVector
{
    const auto& G = standard_partial_molar_gibbs_energies;
    const auto& V = standard_partial_molar_volumes;
    return G - P*V;
}

auto ThermoProperties::standardPartialMolarHeatCapacitiesConstP() const -> ThermoVector
{
    return standard_partial_molar_heat_capacities_cp;
}

auto ThermoProperties::standardPartialMolarHeatCapacitiesConstV() const -> ThermoVector
{
    return standard_partial_molar_heat_capacities_cv;
}


SpeciesThermoProperties::SpeciesThermoProperties()
{}

auto SpeciesThermoProperties::temperature() const -> double
{
    return T;
}

auto SpeciesThermoProperties::pressure() const -> double
{
    return P;
}

auto SpeciesThermoProperties::standardPartialMolarGibbsEnergy() const -> ThermoScalar
{
    return standard_partial_molar_gibbs_energy;
}

auto SpeciesThermoProperties::standardPartialMolarEnthalpy() const -> ThermoScalar
{
    return standard_partial_molar_enthalpy;
}

auto SpeciesThermoProperties::standardPartialMolarVolume() const -> ThermoScalar
{
    return standard_partial_molar_volume;
}

auto SpeciesThermoProperties::standardPartialMolarEntropy() const -> ThermoScalar
{
    const auto& G = standard_partial_molar_gibbs_energy;
    const auto& H = standard_partial_molar_enthalpy;
    return (H - G)/T;

}

auto SpeciesThermoProperties::standardPartialMolarInternalEnergy() const -> ThermoScalar
{
    const auto& H = standard_partial_molar_enthalpy;
    const auto& V = standard_partial_molar_volume;
    return H - P*V;
}

auto SpeciesThermoProperties::standardPartialMolarHelmholtzEnergy() const -> ThermoScalar
{
    const auto& G = standard_partial_molar_gibbs_energy;
    const auto& V = standard_partial_molar_volume;
    return G - P*V;
}

auto SpeciesThermoProperties::standardPartialMolarHeatCapacityConstP() const -> ThermoScalar
{
    return standard_partial_molar_heat_capacity_cp;
}

auto SpeciesThermoProperties::standardPartialMolarHeatCapacityConstV() const -> ThermoScalar
{
    return standard_partial_molar_heat_capacity_cv;
}

} // namespace Reaktoro
