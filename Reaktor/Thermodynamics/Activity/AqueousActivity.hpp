// Reaktor is a C++ library for computational reaction modelling.
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

// Reaktor includes
#include <Reaktor/Thermodynamics/Activity/ActivityUtils.hpp>

namespace Reaktor {

/// A type used to define the function signature of an aqueous activity function
/// @param params An instance of AqueousSolutionState containing the necessary parameters for the activity calculation
/// @return An instance of ChemicalScalar containing the calculated activity and its molar derivatives
/// @see AqueousSolutionState, ChemicalScalar
using AqueousActivity = std::function<ChemicalScalar(const AqueousSolutionState& state)>;

} // namespace Reaktor
