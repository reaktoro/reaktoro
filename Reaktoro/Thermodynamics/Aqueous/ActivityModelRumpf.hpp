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

// Reaktoro includes
#include <Reaktoro/Core/ActivityModel.hpp>

namespace Reaktoro {

/// Return the activity model for a dissolved gas species in an aqueous phase based on Rumpf (1994).
/// **References:**
///  - Rumpf, B., Nicolaisen, H., Ocal, C., & Maurer, G. (1994). Solubility of
///    carbon dioxide in aqueous mixtures of sodium chloride: Experimental
///    results and correlation. Journal of Solution Chemistry, 23(3), 431–448*.
/// @param gas The chemical formula of the dissolved gas.
/// @ingroup ActivityModels
auto ActivityModelRumpf(String gas) -> ActivityModel;

} // namespace Reaktoro
