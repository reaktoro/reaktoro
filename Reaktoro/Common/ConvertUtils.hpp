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

#pragma once

namespace Reaktoro {

/// Enum that defines a temperature unit
enum TemperatureUnit { degC, degK };

/// Enum that defines a pressure unit
enum PressureUnit { Pa, KPa, MPa, bar };

/// Converts temperature from a unit to another
template<TemperatureUnit from, TemperatureUnit to>
double convert(double T) { return T; }

/// Converts pressure from a unit to another
template<PressureUnit from, PressureUnit to>
double convert(double P) { return P; }

/// Converts temperature from celsius to kelvin
template<> inline double convert<degC,degK>(double T) { return T + 273.15; }

/// Converts temperature from kelvin to celsius
template<> inline double convert<degK,degC>(double T) { return T - 273.15; }

/// Converts pressure from pascal to kilo pascal
template<> inline double convert<Pa,KPa>(double P) { return P * 1.0e-3; }

/// Converts pressure from pascal to mega pascal
template<> inline double convert<Pa,MPa>(double P) { return P * 1.0e-6; }

/// Converts pressure from pascal to bar
template<> inline double convert<Pa,bar>(double P) { return P * 1.0e-5; }

/// Converts pressure from kilo pascal to pascal
template<> inline double convert<KPa,Pa>(double P)  { return P * 1.0e+3; }

/// Converts pressure from kilo pascal to mega pascal
template<> inline double convert<KPa,MPa>(double P) { return P * 1.0e-3; }

/// Converts pressure from kilo pascal to bar
template<> inline double convert<KPa,bar>(double P) { return P * 1.0e-2; }

/// Converts pressure from mega pascal to pascal
template<> inline double convert<MPa,Pa>(double P)  { return P * 1.0e+6; }

/// Converts pressure from mega pascal to kilo pascal
template<> inline double convert<MPa,KPa>(double P) { return P * 1.0e+3; }

/// Converts pressure from mega pascal to bar
template<> inline double convert<MPa,bar>(double P) { return P * 1.0e+1; }

/// Converts pressure from bar to pascal
template<> inline double convert<bar,Pa>(double P)  { return P * 1.0e+5; }

/// Converts pressure from bar to kilo pascal
template<> inline double convert<bar,KPa>(double P) { return P * 1.0e+2; }

/// Converts pressure from bar to mega pascal
template<> inline double convert<bar,MPa>(double P) { return P * 1.0e-1; }

} // namespace Reaktoro
