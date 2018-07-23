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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/ScalarTypes.hpp>
#include <Reaktoro/Thermodynamics/Common/StateOfMatter.hpp>

namespace Reaktoro {

/// Calculate the density of water using the Haar--Gallagher--Kell (1984) equation of state
/// @param T The temperature of water (in units of K)
/// @param P The pressure of water (in units of Pa)
/// @return The density of liquid water (in units of kg/m3)
auto waterDensityHGK(Temperature T, Pressure P, StateOfMatter stateofmatter) -> ThermoScalar;

/// Calculate the density of water using the Wagner and Pruss (1995) equation of state
/// @param T The temperature of water (in units of K)
/// @param P The pressure of water (in units of Pa)
/// @return The density of liquid water (in units of kg/m3)
auto waterDensityWagnerPruss(Temperature T, Pressure P, StateOfMatter stateofmatter) -> ThermoScalar;

/// Calculate the density of liquid water using the Haar--Gallagher--Kell (1984) equation of state
/// @param T The temperature of water (in units of K)
/// @param P The pressure of water (in units of Pa)
/// @return The density of liquid water (in units of kg/m3)
auto waterLiquidDensityHGK(Temperature T, Pressure P) -> ThermoScalar;

/// Calculate the density of liquid water using the Wagner and Pruss (1995) equation of state
/// @param T The temperature of water (in units of K)
/// @param P The pressure of water (in units of Pa)
/// @return The density of liquid water (in units of kg/m3)
auto waterLiquidDensityWagnerPruss(Temperature T, Pressure P) -> ThermoScalar;

/// Calculate the density of vapor water using the Haar--Gallagher--Kell (1984) equation of state
/// @param T The temperature of water (in units of K)
/// @param P The pressure of water (in units of Pa)
/// @return The density of water (in units of kg/m3)
auto waterVaporDensityHGK(Temperature T, Pressure P) -> ThermoScalar;

/// Calculate the density of vapor water using the Wagner and Pruss (1995) equation of state
/// @param T The temperature of water (in units of K)
/// @param P The pressure of water (in units of Pa)
/// @return The density of water (in units of kg/m3)
auto waterVaporDensityWagnerPruss(Temperature T, Pressure P) -> ThermoScalar;

/// Calculate the pressure of water using the Haar--Gallagher--Kell (1984) equation of state
/// @param T The temperature of water (in units of K)
/// @param D The density of water (in units of kg/m3)
/// @return The pressure of water (in units of Pa)
auto waterPressureHGK(Temperature T, ThermoScalar D) -> ThermoScalar;

/// Calculate the pressure of water using the Wagner and Pruss (1995) equation of state
/// @param T The temperature of water (in units of K)
/// @param D The density of water (in units of kg/m3)
/// @return The pressure of water (in units of Pa)
auto waterPressureWagnerPruss(Temperature T, ThermoScalar D) -> ThermoScalar;

/// Calculate the saturated pressure of water using the Wagner and Pruss (1995) equation of state
/// @param T The temperature of water (in units of K)
/// @return The saturated pressure of water (in units of Pa)
auto waterSaturatedPressureWagnerPruss(Temperature T) -> ThermoScalar;

/// Calculate the saturated liquid-density of water using the Wagner and Pruss (1995) equation of state
/// @param T The temperature of water (in units of K)
/// @return The saturated liquid-density of water (in units of kg/m3)
auto waterSaturatedLiquidDensityWagnerPruss(Temperature T) -> ThermoScalar;

/// Calculate the saturated vapour-density of water using the Wagner and Pruss (1995) equation of state
/// @param T The temperature of water (in units of K)
/// @return The saturated vapour-density of water (in units of kg/m3)
auto waterSaturatedVapourDensityWagnerPruss(Temperature T) -> ThermoScalar;

} // namespace Reaktoro
