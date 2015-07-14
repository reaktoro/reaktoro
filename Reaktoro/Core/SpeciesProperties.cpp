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

#include "SpeciesProperties.hpp"

namespace Reaktoro {

SpeciesProperties::SpeciesProperties()
{

}

auto SpeciesProperties::temperature() const -> double
{
    return T.val;
}

auto SpeciesProperties::pressure() const -> double
{
    return P.val;
}

auto SpeciesProperties::standardPartialMolarGibbsEnergy() const -> ThermoScalar
{
    return standard_partial_molar_gibbs_energy;
}

auto SpeciesProperties::standardPartialMolarEnthalpy() const -> ThermoScalar
{
    return standard_partial_molar_enthalpy;
}

auto SpeciesProperties::standardPartialMolarVolume() const -> ThermoScalar
{
    return standard_partial_molar_volume;
}

auto SpeciesProperties::standardPartialMolarEntropy() const -> ThermoScalar
{
    const auto& G = standard_partial_molar_gibbs_energy;
    const auto& H = standard_partial_molar_enthalpy;
    return (H - G)/T;

}

auto SpeciesProperties::standardPartialMolarInternalEnergy() const -> ThermoScalar
{
    const auto& H = standard_partial_molar_enthalpy;
    const auto& V = standard_partial_molar_volume;
    return H - P*V;
}

auto SpeciesProperties::standardPartialMolarHelmholtzEnergy() const -> ThermoScalar
{
    const auto& G = standard_partial_molar_gibbs_energy;
    const auto& V = standard_partial_molar_volume;
    return G - P*V;
}

auto SpeciesProperties::standardPartialMolarHeatCapacityConstP() const -> ThermoScalar
{
    return standard_partial_molar_heat_capacity_cp;
}

auto SpeciesProperties::standardPartialMolarHeatCapacityConstV() const -> ThermoScalar
{
    return standard_partial_molar_heat_capacity_cv;
}

} // namespace Reaktoro
