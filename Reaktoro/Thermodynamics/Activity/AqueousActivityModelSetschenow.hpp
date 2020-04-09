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
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityModel.hpp>

namespace Reaktoro {

/// Create the aqueous activity function of a neutral species based on the Setschenow model.
/// @param mixture The aqueous mixture instance containing the aqueous species
/// @param b The Setschenow constant
/// @return The aqueous activity function of the aqueous species
/// @see AqueousMixture, AqueousActivityModel
auto aqueousActivityModelSetschenow(const AqueousMixture& mixture, double b) -> AqueousActivityModel;

} // namespace Reaktoro
