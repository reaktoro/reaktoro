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
#include <Reaktor/Thermodynamics/Activity/AqueousActivity.hpp>

namespace Reaktor {

/// Create the aqueous activity function of a neutral species based on the Setschenow model
///
/// @param species The name of the aqueous neutral species
/// @param solution The aqueous solution instance containing the aqueous species
/// @param b The Setschenow constant
/// @return The aqueous activity function of the aqueous species
/// @see AqueousSolution, AqueousActivity
auto aqueousActivitySetschenow(const std::string& species, const AqueousSolution& solution, double b) -> AqueousActivity;

} // namespace Reaktor
