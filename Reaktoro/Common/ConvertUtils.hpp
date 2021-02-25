// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

/// Converts temperature from celsius to kelvin
template<typename Scalar> auto convertCelsiusToKelvin(Scalar T) -> decltype(T + 273.15) { return T + 273.15; }

/// Converts temperature from kelvin to celsius
template<typename Scalar> auto convertKelvinToCelsius(Scalar T) -> decltype(T - 273.15) { return T - 273.15; }

/// Converts pressure from pascal to kilo pascal
template<typename Scalar> auto convertPascalToKiloPascal(Scalar P) -> decltype(P * 1.0e-3) { return P * 1.0e-3; }

/// Converts pressure from pascal to mega pascal
template<typename Scalar> auto convertPascalToMegaPascal(Scalar P) -> decltype(P * 1.0e-6) { return P * 1.0e-6; }

/// Converts pressure from pascal to bar
template<typename Scalar> auto convertPascalToBar(Scalar P) -> decltype(P * 1.0e-5) { return P * 1.0e-5; }

/// Converts pressure from kilo pascal to pascal
template<typename Scalar> auto convertKiloPascalToPascal(Scalar P)  -> decltype(P * 1.0e+3) { return P * 1.0e+3; }

/// Converts pressure from kilo pascal to mega pascal
template<typename Scalar> auto convertKiloPascalToMegaPascal(Scalar P) -> decltype(P * 1.0e-3) { return P * 1.0e-3; }

/// Converts pressure from kilo pascal to bar
template<typename Scalar> auto convertKiloPascalToBar(Scalar P) -> decltype(P * 1.0e-2) { return P * 1.0e-2; }

/// Converts pressure from mega pascal to pascal
template<typename Scalar> auto convertMegaPascalToPascal(Scalar P)  -> decltype(P * 1.0e+6) { return P * 1.0e+6; }

/// Converts pressure from mega pascal to kilo pascal
template<typename Scalar> auto convertMegaPascalToKiloPascal(Scalar P) -> decltype(P * 1.0e+3) { return P * 1.0e+3; }

/// Converts pressure from mega pascal to bar
template<typename Scalar> auto convertMegaPascalToBar(Scalar P) -> decltype(P * 1.0e+1) { return P * 1.0e+1; }

/// Converts pressure from bar to pascal
template<typename Scalar> auto convertBarToPascal(Scalar P)  -> decltype(P * 1.0e+5) { return P * 1.0e+5; }

/// Converts pressure from bar to kilo pascal
template<typename Scalar> auto convertBarToKiloPascal(Scalar P) -> decltype(P * 1.0e+2) { return P * 1.0e+2; }

/// Converts pressure from bar to mega pascal
template<typename Scalar> auto convertBarToMegaPascal(Scalar P) -> decltype(P * 1.0e-1) { return P * 1.0e-1; }

/// Convert pressure from bar to atm
template<typename Scalar> auto convertBarToAtm(Scalar P) -> decltype(P * 0.986923267) { return P * 0.986923267; }

/// Converts volume from cm3 to m3
template<typename Scalar> auto convertCubicCentimeterToCubicMeter(Scalar V) -> decltype(V * 1.0e-6) { return V * 1.0e-6; }

/// Converts volume from m3 to cm3
template<typename Scalar> auto convertCubicMeterToCubicCentimeter(Scalar V) -> decltype(V * 1.0e+6) { return V * 1.0e+6; }

/// Converts volume from m3 to liter
template<typename Scalar> auto convertCubicMeterToLiter(Scalar V) -> decltype(V * 1.0e+3) { return V * 1.0e+3; }

/// Converts volume from liter to m3
template<typename Scalar> auto convertLiterToCubicMeter(Scalar V) -> decltype(V * 1.0e-3) { return V * 1.0e-3; }

} // namespace Reaktoro
