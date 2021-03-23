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

// C++ includes
#include <limits>

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

/// The value of ln(10)
constexpr auto ln10 = 2.30258509299404590109361379290930926799774169921875;

/// The value of infinity
constexpr auto inf = std::numeric_limits<double>::infinity();

/// The value of NaN
constexpr auto NaN = std::numeric_limits<double>::quiet_NaN();

/// The value of machine precision epsilon
constexpr auto epsilon = std::numeric_limits<double>::epsilon();

} // namespace Reaktoro
