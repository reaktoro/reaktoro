// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Core/StateOfMatter.hpp>

namespace Reaktoro {

/// Calculate the density of water using the Haar--Gallagher--Kell (1984) equation of state
/// @param T The temperature of water (in K)
/// @param P The pressure of water (in Pa)
/// @param stateofmatter The state of matter of water
/// @return The density of liquid water (in kg/m3)
auto waterDensityHGK(real const& T, real const& P, StateOfMatter stateofmatter) -> real;

/// Calculate the density of water using the Wagner and Pruss (1995) equation of state
/// @param T The temperature of water (in K)
/// @param P The pressure of water (in Pa)
/// @param stateofmatter The state of matter of water
/// @return The density of liquid water (in kg/m3)
auto waterDensityWagnerPruss(real const& T, real const& P, StateOfMatter stateofmatter) -> real;

/// Calculate the density of liquid water using the Haar--Gallagher--Kell (1984) equation of state
/// @param T The temperature of water (in K)
/// @param P The pressure of water (in Pa)
/// @return The density of liquid water (in kg/m3)
auto waterLiquidDensityHGK(real const& T, real const& P) -> real;

/// Calculate the density of liquid water using the Wagner and Pruss (1995) equation of state
/// @param T The temperature of water (in K)
/// @param P The pressure of water (in Pa)
/// @return The density of liquid water (in kg/m3)
auto waterLiquidDensityWagnerPruss(real const& T, real const& P) -> real;

/// Calculate the density of vapor water using the Haar--Gallagher--Kell (1984) equation of state
/// @param T The temperature of water (in K)
/// @param P The pressure of water (in Pa)
/// @return The density of water (in kg/m3)
auto waterVaporDensityHGK(real const& T, real const& P) -> real;

/// Calculate the density of vapor water using the Wagner and Pruss (1995) equation of state
/// @param T The temperature of water (in K)
/// @param P The pressure of water (in Pa)
/// @return The density of water (in kg/m3)
auto waterVaporDensityWagnerPruss(real const& T, real const& P) -> real;

/// Calculate the pressure of water using the Haar--Gallagher--Kell (1984) equation of state
/// @param T The temperature of water (in K)
/// @param D The density of water (in kg/m3)
/// @return The pressure of water (in Pa)
auto waterPressureHGK(real const& T, real const& D) -> real;

/// Calculate the pressure of water using the Wagner and Pruss (1995) equation of state
/// @param T The temperature of water (in K)
/// @param D The density of water (in kg/m3)
/// @return The pressure of water (in Pa)
auto waterPressureWagnerPruss(real const& T, real const& D) -> real;

/// Calculate the saturation pressure of water using the Wagner and Pruss (1995) equation of state
/// @param T The temperature of water (in K)
/// @return The saturation pressure of water (in Pa)
auto waterSaturationPressureWagnerPruss(real const& T) -> real;

/// Calculate the saturation liquid-density of water using the Wagner and Pruss (1995) equation of state
/// @param T The temperature of water (in K)
/// @return The saturation liquid-density of water (in kg/m3)
auto waterSaturationLiquidDensityWagnerPruss(real const& T) -> real;

/// Calculate the saturation vapour-density of water using the Wagner and Pruss (1995) equation of state
/// @param T The temperature of water (in K)
/// @return The saturation vapour-density of water (in kg/m3)
auto waterSaturationVapourDensityWagnerPruss(real const& T) -> real;

/// DEPRECATED (use @ref waterSaturationPressureWagnerPruss)
auto waterSaturatedPressureWagnerPruss(real const& T) -> real;

/// DEPRECATED (use @ref waterSaturationLiquidDensityWagnerPruss)
auto waterSaturatedLiquidDensityWagnerPruss(real const& T) -> real;

/// DEPRECATED (use @ref waterSaturationVapourDensityWagnerPruss)
auto waterSaturatedVapourDensityWagnerPruss(real const& T) -> real;

} // namespace Reaktoro
