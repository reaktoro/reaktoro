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

// Reaktoro includes
#include <Reaktoro/Common/Real.hpp>

namespace Reaktoro {

// Forward declarations
struct WaterElectroState;
struct WaterThermoState;

/// Calculate the electrostatic state of water using the model of Johnson and Norton (1991).
/// Reference:
/// - Johnson, J. W., Norton, D. (1991). Critical phenomena in hydrothermal systems; state,
///   thermodynamic, electrostatic, and transport properties of H2O in the critical region.
///   American Journal of Science, 291(6), 541–648. [doi](http://doi.org/10.2475/ajs.291.6.541)
auto waterElectroStateJohnsonNorton(real T, real P, const WaterThermoState& wts) -> WaterElectroState;

} // namespace Reaktoro
