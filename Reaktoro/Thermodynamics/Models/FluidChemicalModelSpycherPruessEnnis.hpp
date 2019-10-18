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
class FluidMixture;

/// Return a chemical model function for a gaseous phase based on the Spycher et al. (2003) model.
/// This model only supports a gaseous phase with species `CO2(g)` and `H2O(g)`.
/// The model is documented in: *Spycher, N., Pruess, K., Ennis-King, J. (2003). CO2-H2O mixtures in the
/// geological sequestration of CO2. I. Assessment and calculation of mutual solubilities from 12 to 100°C
/// and up to 600 bar. Geochimica et Cosmochimica Acta, 67(16), 3015–3031*.
/// @param mixture The gaseous mixture instance
/// @see FluidMixture, PhaseChemicalModel
auto fluidChemicalModelSpycherPruessEnnis(const FluidMixture& mixture)->PhaseChemicalModel;

} // namespace Reaktoro
