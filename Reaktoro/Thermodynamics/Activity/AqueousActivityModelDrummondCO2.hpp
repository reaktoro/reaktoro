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

/// Create the function for the calculation of ln activity coefficient of CO<sub>2</sub>(aq), based on the model of Drummond (1981).
/// The model is documented in the paper Drummond, S. E. (1981). Boiling and mixing of hydrothermal fluids:
/// chemical effects on mineral precipitation. Pennsylvania State University.
/// @param mixture The aqueous mixture instance
/// @return The aqueous activity function of species `CO2(aq)`
/// @see AqueousMixture, AqueousActivityModel
auto aqueousActivityModelDrummondCO2(const AqueousMixture& mixture) -> AqueousActivityModel;

} // namespace Reaktoro
