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

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

/// Parse a reaction equation
auto parseReaction(const String& equation) -> Pairs<String, double>;

/// Parse a formatted string containing pairs of numbers and strings.
/// See below several examples of parsing different formatted strings:
/// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// using namespace Reaktoro;
/// auto pairs0 = parseNumberStringPairs("2:H 1:O");
/// auto pairs1 = parseNumberStringPairs("1:Ca 2:Cl");
/// auto pairs2 = parseNumberStringPairs("1:Mg 1:C 3:O");
/// auto pairs3 = parseNumberStringPairs("1:Ca 1:Mg 2:C 6:O");
/// auto pairs4 = parseNumberStringPairs("3:Fe 2:Al 3:Si 12:O");
/// auto pairs5 = parseNumberStringPairs("1:Na+ 1:Cl-");
/// auto pairs6 = parseNumberStringPairs("1:Ca++ 1:CO3--");
/// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
auto parseNumberStringPairs(const String& str) -> Pairs<String, double>;

} // namespace Reaktoro
