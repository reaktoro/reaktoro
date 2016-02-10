// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

// C++ includes
#include <map>
#include <string>
#include <vector>
#include <istream>

namespace Reaktoro {

// Forward declarations
struct Equilibrium;

/// Return all compounds (e.g., species and titrant names) it can find in an Equilibrium instance.
/// This method is used to setup an automatic chemical system composed of only
/// aqueous species.
auto collectCompounds(const Equilibrium& e) -> std::vector<std::string>;

} // namespace Reaktoro
