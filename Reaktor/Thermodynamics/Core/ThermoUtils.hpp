// Reaktor is a C++ library for computational reaction modelling.
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
#include <vector>

// Reaktor includes
#include <Reaktor/Common/ThermoScalar.hpp>

namespace Reaktor {

// Forward declarations
class Database;

auto standardGibbsEnergy(double T, double P, std::string species, const Database& database) -> ThermoScalar;

auto standardHelmholtzEnergy(double T, double P, std::string species, const Database& database) -> ThermoScalar;

auto standardInternalEnergy(double T, double P, std::string species, const Database& database) -> ThermoScalar;

auto standardEnthalpy(double T, double P, std::string species, const Database& database) -> ThermoScalar;

auto standardEntropy(double T, double P, std::string species, const Database& database) -> ThermoScalar;

auto standardVolume(double T, double P, std::string species, const Database& database) -> ThermoScalar;

auto standardHeatCapacityCp(double T, double P, std::string species, const Database& database) -> ThermoScalar;

} // namespace Reaktor
