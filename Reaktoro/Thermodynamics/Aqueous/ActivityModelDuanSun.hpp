// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Core/ActivityModel.hpp>

namespace Reaktoro {

/// Return the activity model for a dissolved gas species in an aqueous phase based on Duan and Sun (2003).
/// **References:**
///  - Duan, Z., Sun, R. (2003). An improved model calculating CO2 solubility
///    in pure water and aqueous NaCl mixtures from 273 to 533 K and from 0 to
///    2000 bar. Chemical Geology, 193(3-4), 257–271
/// Construct a ActivityModelDuanSun object with given dissolved gas formula.
/// @param gas The chemical formula of the dissolved gas.
/// @ingroup Thermodynamics
auto ActivityModelDuanSun(String gas) -> ActivityModelGenerator;

} // namespace Reaktoro
