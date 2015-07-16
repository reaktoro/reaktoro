// Reaktoro is a C++ library for computational reaction modelling.
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

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Models/PhaseChemicalModel.hpp>

namespace Reaktoro {

// Forward declarations
class GaseousMixture;

/// Return a chemical model function for a gaseous phase based on the Spycher et al. (2003) model.
/// This model only supports a gaseous phase with species `CO2(g)` and `H2O(g)`.
/// The model is documented in: *Spycher, N., Pruess, K., Ennis-King, J. (2003). CO2-H2O mixtures in the
/// geological sequestration of CO2. I. Assessment and calculation of mutual solubilities from 12 to 100°C
/// and up to 600 bar. Geochimica et Cosmochimica Acta, 67(16), 3015–3031*.
/// @param mixture The gaseous mixture instance
/// @see GaseousMixture, PhaseChemicalModel
auto gaseousChemicalModelSpycherEtAl2003(const GaseousMixture& mixture) -> PhaseChemicalModel;

} // namespace Reaktoro
