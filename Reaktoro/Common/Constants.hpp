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

/// The universal gas constant (in units of J/(mol*K))
constexpr auto universalGasConstant = 8.3144621;

/// The Faraday constant (in units of C/mol)
constexpr auto faradayConstant = 96485.3329;

/// The constant factor that converts joule to calorie
constexpr auto jouleToCalorie = 0.239005736;

/// The constant factor that converts calorie to joule
constexpr auto calorieToJoule = 4.184;

/// The conversion factor from bar to pascal
constexpr auto barToPascal = 1.0e+05;

/// The conversion factor from cubic centimeters to cubic meters
constexpr auto cubicCentimeterToCubicMeter = 1.0e-06;

/// The conversion factor from cubic meters to cubic centimeters
constexpr auto cubicMeterToCubicCentimeter = 1.0e+06;

} // namespace Reaktoro
