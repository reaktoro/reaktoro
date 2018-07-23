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

namespace Reaktoro {

/// The universal gas constant (in units of J/(mol*K))
const double universalGasConstant = 8.3144621;

/// The Faraday constant (in units of C/mol)
const double faradayConstant = 96485.3329;

/// The constant factor that converts joule to calorie
const double jouleToCalorie = 0.239005736;

/// The constant factor that converts calorie to joule
const double calorieToJoule = 4.184;

/// The conversion factor from bar to pascal
const double barToPascal = 1.0e+05;

/// The conversion factor from cubic centimeters to cubic meters
const double cubicCentimeterToCubicMeter = 1.0e-06;

/// The conversion factor from cubic meters to cubic centimeters
const double cubicMeterToCubicCentimeter = 1.0e+06;

} // namespace Reaktoro
