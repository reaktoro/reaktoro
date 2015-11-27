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
template<typename Scalar> Scalar convertCelsiusToKelvin(Scalar T) { return T + 273.15; }

/// Converts temperature from kelvin to celsius
template<typename Scalar> Scalar convertKelvinToCelsius(Scalar T) { return T - 273.15; }

/// Converts pressure from pascal to kilo pascal
template<typename Scalar> Scalar convertPascalToKiloPascal(Scalar P) { return P * 1.0e-3; }

/// Converts pressure from pascal to mega pascal
template<typename Scalar> Scalar convertPascalToMegaPascal(Scalar P) { return P * 1.0e-6; }

/// Converts pressure from pascal to bar
template<typename Scalar> Scalar convertPascalToBar(Scalar P) { return P * 1.0e-5; }

/// Converts pressure from kilo pascal to pascal
template<typename Scalar> Scalar convertKiloPascalToPascal(Scalar P)  { return P * 1.0e+3; }

/// Converts pressure from kilo pascal to mega pascal
template<typename Scalar> Scalar convertKiloPascalToMegaPascal(Scalar P) { return P * 1.0e-3; }

/// Converts pressure from kilo pascal to bar
template<typename Scalar> Scalar convertKiloPascalToBar(Scalar P) { return P * 1.0e-2; }

/// Converts pressure from mega pascal to pascal
template<typename Scalar> Scalar convertMegaPascalToPascal(Scalar P)  { return P * 1.0e+6; }

/// Converts pressure from mega pascal to kilo pascal
template<typename Scalar> Scalar convertMegaPascalToKiloPascal(Scalar P) { return P * 1.0e+3; }

/// Converts pressure from mega pascal to bar
template<typename Scalar> Scalar convertMegaPascalToBar(Scalar P) { return P * 1.0e+1; }

/// Converts pressure from bar to pascal
template<typename Scalar> Scalar convertBarToPascal(Scalar P)  { return P * 1.0e+5; }

/// Converts pressure from bar to kilo pascal
template<typename Scalar> Scalar convertBarToKiloPascal(Scalar P) { return P * 1.0e+2; }

/// Converts pressure from bar to mega pascal
template<typename Scalar> Scalar convertBarToMegaPascal(Scalar P) { return P * 1.0e-1; }

/// Convert pressure from bar to atm
template<typename Scalar> Scalar convertBarToAtm(Scalar P) { return P * 0.986923267; }

/// Converts volume from cm3 to m3
template<typename Scalar> Scalar convertCubicCentimeterToCubicMeter(Scalar V) { return V * 1.0e-6; }

/// Converts volume from m3 to cm3
template<typename Scalar> Scalar convertCubicMeterToCubicCentimeter(Scalar V) { return V * 1.0e+6; }

/// Converts volume from m3 to liter
template<typename Scalar> Scalar convertCubicMeterToLiter(Scalar V) { return V * 1.0e+3; }

/// Converts volume from liter to m3
template<typename Scalar> Scalar convertLiterToCubicMeter(Scalar V) { return V * 1.0e-3; }

} // namespace Reaktoro
