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

// C++ includes
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Common/Index.hpp>

namespace Reaktor {

/// A type used to define a chemical element and its attributes
class Element
{
public:
	/// Construct a default Element instance
	Element();

	/// Get the name of the element
	auto name() const -> std::string;
};

/// A type used to define a list of Element instances
typedef std::vector<Element> ElementList;

/// Compare two Element instances
auto operator<(const Element& lhs, const Element& rhs) -> bool;

} // namespace Reaktor
