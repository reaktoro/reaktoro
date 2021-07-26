// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

struct WaterThermoProps
{
    /// The temperature of water (in units of K)
    real temperature = {};

    /// The specific volume of water (in units of m3/kg)
    real volume = {};

    /// The specific entropy of water (in units of J/(kg*K))
    real entropy = {};

    /// The specific Helmholtz free energy of water (in units of J/kg)
    real helmholtz = {};

    /// The specific internal energy of water (in units of J/kg)
    real internal_energy = {};

    /// The specific enthalpy of water (in units of J/kg)
    real enthalpy = {};

    /// The specific Gibbs free energy of water (in units of J/kg)
    real gibbs = {};

    /// The specific isochoric heat capacity of water (in units of J/(kg*K))
    real cv = {};

    /// The specific isobaric heat capacity of water (in units of J/(kg*K))
    real cp = {};

    /// The specific density of water (in units of kg/m3)
    real density = {};

    /// The first-order partial derivative of density with respect to temperature (in units of (kg/m3)/K)
    real densityT = {};

    /// The first-order partial derivative of density with respect to pressure (in units of (kg/m3)/Pa)
    real densityP = {};

    /// The second-order partial derivative of density with respect to temperature (in units of (kg/m3)/(K*K))
    real densityTT = {};

    /// The second-order partial derivative of density with respect to temperature and pressure (in units of (kg/m3)/(K*Pa))
    real densityTP = {};

    /// The second-order partial derivative of density with respect to pressure (in units of (kg/m3)/(Pa*Pa))
    real densityPP = {};

    /// The pressure of water (in units of Pa)
    real pressure = {};

    /// The first-order partial derivative of pressure with respect to temperature (in units of Pa/K)
    real pressureT = {};

    /// The first-order partial derivative of pressure with respect to density (in units of Pa/(kg/m3))
    real pressureD = {};

    /// The second-order partial derivative of pressure with respect to temperature (in units of Pa/(K*K))
    real pressureTT = {};

    /// The second-order partial derivative of pressure with respect to temperature and density (in units of Pa/(K*kg/m3))
    real pressureTD = {};

    /// The second-order partial derivative of pressure with respect to density (in units of Pa/((kg/m3)*(kg/m3)))
    real pressureDD = {};
};

} // namespace Reaktoro
