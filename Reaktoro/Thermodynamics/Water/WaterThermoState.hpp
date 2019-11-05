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
#include <Reaktoro/Common/ThermoScalar.hpp>

namespace Reaktoro {

struct WaterThermoState
{
    /// The temperature of water (in units of K)
    ThermoScalar temperature;

    /// The specific volume of water (in units of m3/kg)
    ThermoScalar volume;

    /// The specific entropy of water (in units of J/(kg*K))
    ThermoScalar entropy;

    /// The specific Helmholtz free energy of water (in units of J/kg)
    ThermoScalar helmholtz;

    /// The specific internal energy of water (in units of J/kg)
    ThermoScalar internal_energy;

    /// The specific enthalpy of water (in units of J/kg)
    ThermoScalar enthalpy;

    /// The specific Gibbs free energy of water (in units of J/kg)
    ThermoScalar gibbs;

    /// The specific isochoric heat capacity of water (in units of J/(kg*K))
    ThermoScalar cv;

    /// The specific isobaric heat capacity of water (in units of J/(kg*K))
    ThermoScalar cp;

    /// The specific density of water (in units of kg/m3)
    ThermoScalar density;

    /// The first-order partial derivative of density with respect to temperature (in units of (kg/m3)/K)
    ThermoScalar densityT;

    /// The first-order partial derivative of density with respect to pressure (in units of (kg/m3)/Pa)
    ThermoScalar densityP;

    /// The second-order partial derivative of density with respect to temperature (in units of (kg/m3)/(K*K))
    ThermoScalar densityTT;

    /// The second-order partial derivative of density with respect to temperature and pressure (in units of (kg/m3)/(K*Pa))
    ThermoScalar densityTP;

    /// The second-order partial derivative of density with respect to pressure (in units of (kg/m3)/(Pa*Pa))
    ThermoScalar densityPP;

    /// The pressure of water (in units of Pa)
    ThermoScalar pressure;

    /// The first-order partial derivative of pressure with respect to temperature (in units of Pa/K)
    ThermoScalar pressureT;

    /// The first-order partial derivative of pressure with respect to density (in units of Pa/(kg/m3))
    ThermoScalar pressureD;

    /// The second-order partial derivative of pressure with respect to temperature (in units of Pa/(K*K))
    ThermoScalar pressureTT;

    /// The second-order partial derivative of pressure with respect to temperature and density (in units of Pa/(K*kg/m3))
    ThermoScalar pressureTD;

    /// The second-order partial derivative of pressure with respect to density (in units of Pa/((kg/m3)*(kg/m3)))
    ThermoScalar pressureDD;
};

} // namespace Reaktoro
