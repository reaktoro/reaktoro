// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

// C++ includes
#include <string>

namespace units {

/// Return the slope factor in the linear function that converts a numeric value from a unit to another.
/// @param from The string representing the unit from which the conversion is made
/// @param to The string representing the unit to which the conversion is made
/// @return The slope factor
auto slope(const std::string& from, const std::string& to) -> double;

/// Return the intercept term in the linear function that converts a numeric value from a unit to another.
/// @param from The string representing the unit from which the conversion is made
/// @param to The string representing the unit to which the conversion is made
/// @return The intercept term
auto intercept(const std::string& from, const std::string& to) -> double;

/// Check if two units are convertible among each other
/// @return True if they are convertible, false otherwise
auto convertible(const std::string& from, const std::string& to) -> bool;

/// Convert a numeric value from a unit to another
/// @param value The value
/// @param from The string representing the unit from which the conversion is made
/// @param to The string representing the unit to which the conversion is made
/// @return The converted value
template<typename T>
auto convert(const T& value, const std::string& from, const std::string& to) -> T
{
    return (from == to) ? value : value * slope(from, to) + intercept(from, to);
}

} /* namespace units */
