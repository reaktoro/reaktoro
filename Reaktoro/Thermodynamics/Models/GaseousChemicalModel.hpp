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
#include <Reaktoro/Thermodynamics/Mixtures/GaseousMixture.hpp>

namespace Reaktoro {

/// The signature of an equation of state function that calculates the thermodynamic properties of gaseous phases.
/// @see GaseousMixtureState, ChemicalModel, PhaseChemicalModelResult
using GaseousChemicalModel = std::function<PhaseChemicalModelResult(const GaseousMixtureState&)>;

} // namespace Reaktoro
