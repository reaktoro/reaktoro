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

// Reaktoro includes
#include <Reaktoro/Core/ActivityModel.hpp>

namespace Reaktoro {

/// Return the activity model for gaseous phases formulated in Spycher et al. (2003).
/// This is an activity model for a gaseous phase supporting only the gases
/// H<sub>2</sub>2O(g) and CO<sub>2</sub>(g).
///
/// **References:**
/// - Spycher, N., Pruess, K., Ennis-King, J. (2003). CO2-H2O mixtures in the
///   geological sequestration of CO2. I. Assessment and calculation of mutual
///   solubilities from 12 to 100C and up to 600 bar. Geochimica et
///   Cosmochimica Acta, 67(16), 3015-3031.
auto ActivityModelSpycherPruessEnnis() -> ActivityModel;

} // namespace Reaktoro
