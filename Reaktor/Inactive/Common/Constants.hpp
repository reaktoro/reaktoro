/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

namespace Reaktor {

/// The universal gas constant (in units of J/(mol*K))
const double universalGasConstant = 8.3144621;

/// The constant factor that converts joule to calorie
const double jouleToCalorie = 0.239005736;

/// The constant factor that converts calorie to joule
const double calorieToJoule = 4.184;

/// The conversion factor from bar to pascal
const double barToPascal = 1.0e-05;

/// The conversion factor from cubic centimeters to cubic meters
const double cubicCentimeterToCubicMeter = 1.0e-06;

/// The conversion factor from cubic meters to cubic centimeters
const double cubicMeterToCubicCentimeter = 1.0e+06;

} /* namespace Reaktor */
