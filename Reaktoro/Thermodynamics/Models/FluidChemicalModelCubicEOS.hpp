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

/// Set the chemical model of the phase with the van der Waals equation of state.
/// Reference: *van der Waals, J.D. (1910). Nobel Lectures in Physics. pp. 254-265*.
auto fluidChemicalModelVanDerWaals(const FluidMixture& mixture)->PhaseChemicalModel;

/// Set the chemical model of the phase with the Redlich-Kwong equation of state.
/// Reference: *Redlich, O., Kwong, J.N.S. (1949). On The Thermodynamics of Solutions. Chem. Rev. 44(1) 233–244*.
auto fluidChemicalModelRedlichKwong(const FluidMixture& mixture)->PhaseChemicalModel;

/// Set the chemical model of the phase with the Soave-Redlich-Kwong equation of state.
/// Reference: *Soave, G. (1972). Equilibrium constants from a modified Redlich-Kwong equation of state, Chem. Eng. Sci., 27, 1197-1203*.
auto fluidChemicalModelSoaveRedlichKwong(const FluidMixture& mixture)->PhaseChemicalModel;

/// Set the chemical model of the phase with the Peng-Robinson equation of state.
/// Reference: *Peng, D.Y., Robinson, D.B. (1976). A New Two-Constant Equation of State. Industrial and Engineering Chemistry: Fundamentals 15: 59–64*.
auto fluidChemicalModelPengRobinson(const FluidMixture& mixture)->PhaseChemicalModel;

} // namespace Reaktoro
