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

/// Return a chemical model function for a gaseous phase based on the Spycher and Reed (1988) model.
/// This model only supports the gaseous species `H2O(g)`, `CO2(g)`, and `CH4(g)`.
///
/// Reference: *Spycher, N., Reed, M. (1988). Fugacity coefficients of H2, CO2,
/// CH4, H2O and of H2O--CO2--CH4 mixtures: A virial equation treatment for
/// moderate pressures and temperatures applicable to calculations of
/// hydrothermal boiling. Geochimica et Cosmochimica Acta, 52(3), 739–749*.
/// @param mixture The gaseous mixture instance
/// @see GaseousMixture, GaseousActivityFunction
auto gaseousChemicalModelSpycherReedH2OCO2CH4(const GaseousMixture& mixture) -> PhaseChemicalModel;

} // namespace Reaktoro
