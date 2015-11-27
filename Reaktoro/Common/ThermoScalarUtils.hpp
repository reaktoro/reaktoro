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

// C++ includes
#include <string>

namespace Reaktoro {

// Forward declarations
class Pressure;
class Temperature;

/// Return a ThermoScalar instance for a given temperature value
auto temperature(double val) -> Temperature;

/// Return a ThermoScalar instance for a given temperature value with units
auto temperature(double val, std::string units) -> Temperature;

/// Return a ThermoScalar instance for a given pressure value
auto pressure(double val) -> Pressure;

/// Return a ThermoScalar instance for a given pressure value with units
auto pressure(double val, std::string units) -> Pressure;

} // namespace Reaktoro
