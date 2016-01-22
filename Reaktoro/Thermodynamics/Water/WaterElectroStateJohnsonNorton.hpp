// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

// Reaktoro includes
#include <Reaktoro/Common/ScalarTypes.hpp>

namespace Reaktoro {

// Forward declarations
struct WaterElectroState;
struct WaterThermoState;

// Calculate the electrostatic state of water using the model of Johnson and Norton (1991)
auto waterElectroStateJohnsonNorton(Temperature T, Pressure P, const WaterThermoState& wts) -> WaterElectroState;

} // namespace Reaktoro
