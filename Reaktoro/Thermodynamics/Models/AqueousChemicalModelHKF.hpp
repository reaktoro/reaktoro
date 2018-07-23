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
#include <Reaktoro/Thermodynamics/Models/PhaseChemicalModel.hpp>

namespace Reaktoro {

// Forward declarations
class AqueousMixture;

/// Return an equation of state for an aqueous phase based on HKF model.
/// The HKF model implemented here was documented in:
///   - Helgeson, H. C., Kirkham, D. H., Flowers, G. C. (1981). Theoretical prediction of the thermodynamic behavior
///     of aqueous electrolytes at high pressures and temperatures: IV. Calculation of activity coefficients, osmotic
///     coefficients, and apparent molal and standard and relative partial molal properties to 600°C.
///     American Journal of Science, 281(10), 1249–1516.
/// @param mixture The aqueous mixture
/// @return The equation of state function for the aqueous phase
/// @see AqueousMixture, PhaseChemicalModel
auto aqueousChemicalModelHKF(const AqueousMixture& mixture) -> PhaseChemicalModel;

} // namespace Reaktoro
