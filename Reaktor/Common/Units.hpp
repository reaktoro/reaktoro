/*
 * Units.hpp
 *
 *  Created on: 28 Oct 2013
 *      Author: allan
 */

#pragma once

// C++ includes
#include <string>

namespace units {

/// Convert a numeric value from a unit to another
/// @param value The value
/// @param from The string representing the unit from which the conversion is made
/// @param to The string representing the unit to which the conversion is made
/// @return The converted value
auto convert(double value, const std::string& from, const std::string& to) -> double;

/// Check if two units are convertible among each other
/// @return True if they are convertible, false otherwise
auto convertible(const std::string& from, const std::string& to) -> bool;

} /* namespace units */
