// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

namespace Reaktoro {

/// The molar mass of water in units of kg/mol
const double waterMolarMass = 0.018015268;

/// The critical temperature of water in units of K
const double waterCriticalTemperature = 647.096;

/// The critical pressure of water in units of Pa
const double waterCriticalPressure = 22.064e+06;

/// The critical density of water in units of kg/m3
const double waterCriticalDensity = 322.0;

/// The triple point temperature of water in units of K
const double waterTriplePointTemperature = 273.16;

/// The triple point pressure of water in units of Pa
const double waterTriplePointPressure = 611.655;

/// The triple point liquid-density of water in units of kg/m3
const double waterTriplePointDensityLiquid = 999.793;

/// The triple point vapour-density of water in units of kg/m3
const double waterTriplePointDensityVapour = 0.00485458;

} // namespace Reaktoro
