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
    return T.val;
}

auto ThermoProperties::pressure() const -> double
{
    return P.val;
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

} // namespace Reaktoro
