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

// Forward declarations
class AqueousMixture;

/// Return an equation of state for an aqueous phase based on a Pitzer model.
/// The implementation of this Pitzer model was taken from the following references:
///   1. Harvie, C.E., Møller, N., Weare, J.H. (1984). The prediction of mineral solubilities in natural waters:
///      The Na-K-Mg-Ca-H-Cl-SO4-OH-HCO3-CO3-CO2-H2O system to high ionic strengths at 25°C.
///      Geochimica et Cosmochimica Acta, 48(4), 723–751.
///   2. Harvie, C.E., Weare, J.H. (1980). The prediction of mineral soluhilities in natural waters: the system
///      from zero to high concentration at 25°C. Geochimica et Cosmochimica Acta, 44(7), 981–997.
///   3. Pitzer, K. S. (1975). Thermodynamics of electrolytes. V. effects of higher-order electrostatic terms.
///      Journal of Solution Chemistry, 4(3), 249–265.
/// @param mixture The aqueous mixture
/// @return The equation of state function for the aqueous phase
/// @see AqueousMixture, ActivityPropsFn
auto aqueousChemicalModelPitzerHMW(const AqueousMixture& mixture)-> ActivityPropsFn;

} // namespace Reaktoro
