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

namespace Reaktoro {

/// Used to store thermodynamic properties of water.
struct WaterThermoProps
{
    /// The temperature of water (in K).
    real T;

    /// The specific volume of water (in m3/kg).
    real V;

    /// The specific entropy of water (in J/(kg*K)).
    real S;

    /// The specific Helmholtz free energy of water (in J/kg).
    real A;

    /// The specific internal energy of water (in J/kg).
    real U;

    /// The specific enthalpy of water (in J/kg).
    real H;

    /// The specific Gibbs free energy of water (in J/kg).
    real G;

    /// The specific isochoric heat capacity of water (in J/(kg*K)).
    real Cv;

    /// The specific isobaric heat capacity of water (in J/(kg*K)).
    real Cp;

    /// The specific density of water (in kg/m3).
    real D;

    /// The first-order partial derivative of density with respect to temperature (in (kg/m3)/K).
    real DT;

    /// The first-order partial derivative of density with respect to pressure (in (kg/m3)/Pa).
    real DP;

    /// The second-order partial derivative of density with respect to temperature (in (kg/m3)/(K*K)).
    real DTT;

    /// The second-order partial derivative of density with respect to temperature and pressure (in (kg/m3)/(K*Pa)).
    real DTP;

    /// The second-order partial derivative of density with respect to pressure (in (kg/m3)/(Pa*Pa)).
    real DPP;

    /// The pressure of water (in Pa).
    real P;

    /// The first-order partial derivative of pressure with respect to temperature (in Pa/K).
    real PT;

    /// The first-order partial derivative of pressure with respect to density (in Pa/(kg/m3)).
    real PD;

    /// The second-order partial derivative of pressure with respect to temperature (in Pa/(K*K)).
    real PTT;

    /// The second-order partial derivative of pressure with respect to temperature and density (in Pa/(K*kg/m3)).
    real PTD;

    /// The second-order partial derivative of pressure with respect to density (in Pa/((kg/m3)*(kg/m3))).
    real PDD;

    /// Add another WaterThermoProps object to this.
    auto operator+=(WaterThermoProps const& x) -> WaterThermoProps&;

    /// Subtract another WaterThermoProps object to this.
    auto operator-=(WaterThermoProps const& x) -> WaterThermoProps&;

    /// Multiply this WaterThermoProps object by a scalar.
    auto operator*=(double const& x) -> WaterThermoProps&;

    /// Multiply this WaterThermoProps object by a scalar.
    auto operator*=(real const& x) -> WaterThermoProps&;

    /// Divide this WaterThermoProps object by a scalar.
    auto operator/=(double const& x) -> WaterThermoProps&;

    /// Divide this WaterThermoProps object by a scalar.
    auto operator/=(real const& x) -> WaterThermoProps&;
};

auto operator+(WaterThermoProps const& r) -> WaterThermoProps;
auto operator+(WaterThermoProps&& r) -> WaterThermoProps&&;

auto operator-(WaterThermoProps const& r) -> WaterThermoProps;
auto operator-(WaterThermoProps&& r) -> WaterThermoProps&&;

auto operator+(WaterThermoProps const& l, WaterThermoProps const& r) -> WaterThermoProps;
auto operator+(WaterThermoProps&& l, WaterThermoProps const& r) -> WaterThermoProps&&;
auto operator+(WaterThermoProps const& l, WaterThermoProps&& r) -> WaterThermoProps&&;
auto operator+(WaterThermoProps&& l, WaterThermoProps&& r) -> WaterThermoProps&&;

auto operator-(WaterThermoProps const& l, WaterThermoProps const& r) -> WaterThermoProps;
auto operator-(WaterThermoProps&& l, WaterThermoProps const& r) -> WaterThermoProps&&;
auto operator-(WaterThermoProps const& l, WaterThermoProps&& r) -> WaterThermoProps&&;
auto operator-(WaterThermoProps&& l, WaterThermoProps&& r) -> WaterThermoProps&&;

auto operator*(double const& l, WaterThermoProps const& r) -> WaterThermoProps;
auto operator*(double const& l, WaterThermoProps&& r) -> WaterThermoProps&&;
auto operator*(WaterThermoProps const& l, double const& r) -> WaterThermoProps;
auto operator*(WaterThermoProps&& l, double const& r) -> WaterThermoProps&&;

auto operator*(real const& l, WaterThermoProps const& r) -> WaterThermoProps;
auto operator*(real const& l, WaterThermoProps&& r) -> WaterThermoProps&&;
auto operator*(WaterThermoProps const& l, real const& r) -> WaterThermoProps;
auto operator*(WaterThermoProps&& l, real const& r) -> WaterThermoProps&&;

auto operator/(WaterThermoProps const& l, double const& r) -> WaterThermoProps;
auto operator/(WaterThermoProps&& l, double const& r) -> WaterThermoProps&&;

auto operator/(WaterThermoProps const& l, real const& r) -> WaterThermoProps;
auto operator/(WaterThermoProps&& l, real const& r) -> WaterThermoProps&&;

} // namespace Reaktoro
