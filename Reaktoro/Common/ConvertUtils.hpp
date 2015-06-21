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

/// Converts temperature from celsius to kelvin
inline double convertCelsiusToKelvin(double T) { return T + 273.15; }

/// Converts temperature from kelvin to celsius
inline double convertKelvinToCelsius(double T) { return T - 273.15; }

/// Converts pressure from pascal to kilo pascal
inline double convertPascalToKiloPascal(double P) { return P * 1.0e-3; }

/// Converts pressure from pascal to mega pascal
inline double convertPascalToMegaPascal(double P) { return P * 1.0e-6; }

/// Converts pressure from pascal to bar
inline double convertPascalToBar(double P) { return P * 1.0e-5; }

/// Converts pressure from kilo pascal to pascal
inline double convertKiloPascalToPascal(double P)  { return P * 1.0e+3; }

/// Converts pressure from kilo pascal to mega pascal
inline double convertKiloPascalToMegaPascal(double P) { return P * 1.0e-3; }

/// Converts pressure from kilo pascal to bar
inline double convertKiloPascalToBar(double P) { return P * 1.0e-2; }

/// Converts pressure from mega pascal to pascal
inline double convertMegaPascalToPascal(double P)  { return P * 1.0e+6; }

/// Converts pressure from mega pascal to kilo pascal
inline double convertMegaPascalToKiloPascal(double P) { return P * 1.0e+3; }

/// Converts pressure from mega pascal to bar
inline double convertMegaPascalToBar(double P) { return P * 1.0e+1; }

/// Converts pressure from bar to pascal
inline double convertBarToPascal(double P)  { return P * 1.0e+5; }

/// Converts pressure from bar to kilo pascal
inline double convertBarToKiloPascal(double P) { return P * 1.0e+2; }

/// Converts pressure from bar to mega pascal
inline double convertBarToMegaPascal(double P) { return P * 1.0e-1; }

} // namespace Reaktoro
